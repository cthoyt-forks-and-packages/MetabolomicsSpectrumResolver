import math

def fast_cosine_shift(spectrum1,spectrum2,tol,min_match):
    if spectrum1['n_peaks'] == 0 or spectrum2['n_peaks'] == 0:
        return 0.0,[]

    spec1 = spectrum1['normalised_peaks']
    spec2 = spectrum2['normalised_peaks']

    zero_pairs = find_pairs(spec1,spec2,tol,shift=0.0)

    shift = spectrum1['precursor_mz'] - spectrum2['precursor_mz']

    nonzero_pairs = find_pairs(spec1,spec2,tol,shift = shift)

    matching_pairs = zero_pairs + nonzero_pairs

    matching_pairs = sorted(matching_pairs,key = lambda x: x[2], reverse = True)

    used1 = set()
    used2 = set()
    score = 0.0
    used_matches = []
    for m in matching_pairs:
        if not m[0] in used1 and not m[1] in used2:
            score += m[2]
            used1.add(m[0])
            used2.add(m[1])
            used_matches.append(m)
    if len(used_matches) < min_match:
        score = 0.0
    return score,used_matches

def sqrt_normalise(peaks):
    temp = []
    total = 0.0
    for mz,intensity in peaks:
        temp.append((mz,math.sqrt(intensity)))
        total += intensity
    norm_facc = math.sqrt(total)
    normalised_peaks = []
    for mz,intensity in temp:
        normalised_peaks.append((mz,intensity/norm_facc))
    return normalised_peaks

def plot_spectral_alignment(spectrum1,spectrum2,similarity_function,similarity_tolerance,outname,scale = True,**kwargs):
    import pylab as plt
    score,matches = similarity_function(spectrum1,spectrum2,similarity_tolerance,0) # set min_match to zero
    cols = ['r','b','g','k','m','c']
    plt.figure(**kwargs)
    if scale:
        maxi1 = max([intensity for mz,intensity in spectrum1['peaks']])
        maxi2 = max([intensity for mz,intensity in spectrum2['peaks']])
    else:
        maxi1 = 1.0
        maxi2 = 1.0
        
    for mz,intensity in spectrum1['peaks']:
        plt.plot([mz,mz],[0,intensity/maxi1],'k',color=[0.3, 0.3, 0.3])
    for mz,intensity in spectrum2['peaks']:
        plt.plot([mz,mz],[0,-intensity/maxi2],'k',color=[0.3, 0.3, 0.3])
    plt.plot(plt.xlim(),[0,0],'k--',color =[0.6,0.6,0.6])
    colpos = 0
    for pos1,pos2,_ in matches:
        plt.plot([spectrum1['peaks'][pos1][0],spectrum1['peaks'][pos1][0]],[0,spectrum1['peaks'][pos1][1]/maxi1],cols[colpos])
        plt.plot([spectrum2['peaks'][pos2][0],spectrum2['peaks'][pos2][0]],[0,-spectrum2['peaks'][pos2][1]/maxi2],cols[colpos])
        if abs(spectrum1['peaks'][pos1][0] - spectrum2['peaks'][pos2][0]) >= similarity_tolerance:
            plt.plot([spectrum1['peaks'][pos1][0],spectrum2['peaks'][pos2][0]],[spectrum1['peaks'][pos1][1]/maxi1,-spectrum2['peaks'][pos2][1]/maxi2],'k--',color = [0.6,0.6,0.6])
        colpos += 1
        if colpos == len(cols):
            colpos = 0
    plt.title("{:.2f} <-> {:.2f}, Score = {}".format(spectrum1['precursor_mz'],spectrum2['precursor_mz'],score))
    plt.savefig(outname)

def find_pairs(spec1,spec2,tol,shift=0):
    matching_pairs = []
    spec2lowpos = 0
    spec2length = len(spec2)
    
    for idx,(mz,intensity) in enumerate(spec1):
        # do we need to increase the lower idx?
        while spec2lowpos < spec2length and spec2[spec2lowpos][0] + shift < mz - tol:
            spec2lowpos += 1
        if spec2lowpos == spec2length:
            break
        spec2pos = spec2lowpos
        while(spec2pos < spec2length and spec2[spec2pos][0] + shift < mz + tol):
            matching_pairs.append((idx,spec2pos,intensity*spec2[spec2pos][1]))
            spec2pos += 1
        
    return matching_pairs    