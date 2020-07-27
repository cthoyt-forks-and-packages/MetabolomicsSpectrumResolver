import functools
import json
from typing import Tuple

import requests
import spectrum_utils.spectrum as sus

MS2LDA_SERVER = 'http://ms2lda.org/basicviz/'
MOTIFDB_SERVER = 'http://ms2lda.org/motifdb/'
MASSBANK_SERVER = 'https://massbank.us/rest/spectra/'


@functools.lru_cache(100)
def parse_usi_legacy(usi: str) -> Tuple[sus.MsmsSpectrum, str]:
    usi_identifier = usi.lower().split(':')[1]
    if usi_identifier.startswith('gnpstask'):
        return _parse_gnps_task(usi)
    elif usi_identifier.startswith('gnpslibrary'):
        return _parse_gnps_library(usi)
    elif usi_identifier.startswith('ms2ldatask'):
        return _parse_ms2lda(usi)
    elif usi_identifier.startswith('pxd'):
        return _parse_msv_pxd(usi)
    elif usi_identifier.startswith('msv'):
        return _parse_msv_pxd(usi)
    elif usi_identifier.startswith('mtbls'):
        return _parse_mtbls(usi)
    elif usi_identifier.startswith('st'):
        return _parse_metabolomics_workbench(usi)
    elif usi_identifier.startswith('motifdb'):
        return _parse_motifdb(usi)
    elif usi_identifier.startswith('massbank'):
        return _parse_massbank(usi)
    else:
        raise ValueError(f'Unknown USI: {usi}')


# Parse GNPS clustered spectra in Molecular Networking.
def _parse_gnps_task(usi: str) -> Tuple[sus.MsmsSpectrum, str]:
    tokens = usi.split(':')
    task = tokens[1].split('-')[1]
    filename = tokens[2]
    scan = tokens[4]
    request_url = (f'https://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?'
                   f'task={task}&invoke=annotatedSpectrumImageText&block=0&'
                   f'file=FILE->{filename}&scan={scan}&peptide=*..*&'
                   f'force=false&_=1561457932129&format=JSON')
    spectrum_dict = requests.get(request_url).json()
    mz, intensity = zip(*spectrum_dict['peaks'])
    source_link = f'https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task={task}'
    if 'precursor' in spectrum_dict:
        precursor_mz = float(spectrum_dict['precursor'].get('mz', 0))
        charge = int(spectrum_dict['precursor'].get('charge', 1))
    else:
        precursor_mz, charge = 0, 1
    return (sus.MsmsSpectrum(usi, precursor_mz, charge, mz, intensity),
            source_link)


# Parse GNPS library.
def _parse_gnps_library(usi: str) -> Tuple[sus.MsmsSpectrum, str]:
    tokens = usi.split(':')
    identifier = tokens[2]
    request_url = (f'https://gnps.ucsd.edu/ProteoSAFe/SpectrumCommentServlet?'
                   f'SpectrumID={identifier}')
    spectrum_dict = requests.get(request_url).json()
    mz, intensity = zip(*json.loads(
        spectrum_dict['spectruminfo']['peaks_json']))
    source_link = (f'https://gnps.ucsd.edu/ProteoSAFe/'
                   f'gnpslibraryspectrum.jsp?SpectrumID={identifier}')
    return sus.MsmsSpectrum(
        usi, float(spectrum_dict['annotations'][0]['Precursor_MZ']), 1, mz,
        intensity), source_link


# Parse MS2LDA from ms2lda.org.
def _parse_ms2lda(usi: str) -> Tuple[sus.MsmsSpectrum, str]:
    tokens = usi.split(':')
    experiment_id = tokens[1].split('-')[1]
    document_id = tokens[3]
    request_url = (f'{MS2LDA_SERVER}get_doc/?experiment_id={experiment_id}'
                   f'&document_id={document_id}')
    spectrum_dict = json.loads(requests.get(request_url).text)
    mz, intensity = zip(*spectrum_dict['peaks'])
    source_link = f'http://ms2lda.org/basicviz/show_doc/{document_id}/'
    return sus.MsmsSpectrum(usi, float(spectrum_dict['precursor_mz']), 1, mz,
                            intensity), source_link


# Parse MSV or PXD library.
def _parse_msv_pxd(usi: str) -> Tuple[sus.MsmsSpectrum, str]:
    tokens = usi.split(':')
    dataset_identifier = tokens[1]
    filename = tokens[2]
    scan = tokens[4]
    lookup_url = (f'https://massive.ucsd.edu/ProteoSAFe/QuerySpectrum?'
                  f'id=mzspec:{dataset_identifier}:{filename}:scan:{scan}')
    usi_resolvable = False
    for spectrum_file in requests.get(lookup_url).json()['row_data']:
        try:
            usi_resolvable = any(
                spectrum_file['file_descriptor'].lower().endswith(extension)
                for extension in ['mzml', 'mzxml', 'mgf'])
            if usi_resolvable:
                request_url = (f'https://gnps.ucsd.edu/ProteoSAFe/'
                               f'DownloadResultFile?'
                               f'task=4f2ac74ea114401787a7e96e143bb4a1&'
                               f'invoke=annotatedSpectrumImageText&block=0&'
                               f'file=FILE->{spectrum_file["file_descriptor"]}'
                               f'&scan={scan}&peptide=*..*&force=false&'
                               f'format=JSON&uploadfile=True')
                spectrum_dict = requests.get(request_url).json()
                mz, intensity = zip(*spectrum_dict['peaks'])

                precursor_mz = 0
                charge = 1

                try:
                    charge = int(spectrum_dict['precursor']['charge'])
                    precursor_mz = float(spectrum_dict['precursor']['mz'])
                except:
                    pass

                if dataset_identifier.startswith('PXD'):
                    source_link = (
                        f'http://proteomecentral.proteomexchange.org/'
                        f'cgi/GetDataset?ID={dataset_identifier}')
                else:
                    source_link = (f'https://massive.ucsd.edu/ProteoSAFe/'
                                   f'QueryMSV?id={dataset_identifier}')

                return sus.MsmsSpectrum(usi, precursor_mz, charge, mz,
                                        intensity), source_link
        except:
            pass
    if usi_resolvable:
        raise ValueError('Cannot resolve USI')
    raise ValueError('Unsupported/unknown USI')


def _parse_mtbls(usi: str) -> Tuple[sus.MsmsSpectrum, str]:
    tokens = usi.split(':')
    dataset_identifier = tokens[1]
    filename = tokens[2]
    scan = tokens[4]
    for dataset in requests.get('https://massive.ucsd.edu/ProteoSAFe/'
                                'datasets_json.jsp').json()['datasets']:
        if dataset_identifier in dataset['title']:
            source_link = (f'https://www.ebi.ac.uk/'
                           f'metabolights/{dataset_identifier}')
            return _parse_msv_pxd(f'mzspec:{dataset["dataset"]}:{filename}:'
                                  f'scan:{scan}')[0], source_link
    raise ValueError('Unsupported/unknown USI')


def _parse_metabolomics_workbench(usi: str) -> Tuple[sus.MsmsSpectrum, str]:
    tokens = usi.split(':')
    dataset_identifier = tokens[1]
    filename = tokens[2]
    scan = tokens[4]
    for dataset in requests.get('https://massive.ucsd.edu/ProteoSAFe/'
                                'datasets_json.jsp').json()['datasets']:
        if dataset_identifier in dataset['title']:
            source_link = (f'https://www.metabolomicsworkbench.org/'
                           f'data/DRCCMetadata.php?Mode=Study&StudyID=/'
                           f'{dataset_identifier}')
            return _parse_msv_pxd(f'mzspec:{dataset["dataset"]}:{filename}:'
                                  f'scan:{scan}')[0], source_link
    raise ValueError('Unsupported/unknown USI')


# Parse MOTIFDB from ms2lda.org.
def _parse_motifdb(usi: str) -> Tuple[sus.MsmsSpectrum, str]:
    # E.g. mzspec:MOTIFDB:motif:motif_id.
    tokens = usi.split(':')
    motif_id = tokens[3]
    request_url = f'{MOTIFDB_SERVER}get_motif/{motif_id}'
    mz, intensity = zip(*json.loads(requests.get(request_url).text))
    source_link = f'http://ms2lda.org/motifdb/motif/{motif_id}/'
    return sus.MsmsSpectrum(usi, 0, 1, mz, intensity), source_link


# Parse MassBank entry.
def _parse_massbank(usi: str) -> Tuple[sus.MsmsSpectrum, str]:
    # E.g. mzspec:MASSBANK:motif:motif_id.
    tokens = usi.split(':')
    massbank_id = tokens[2]
    request_url = f'{MASSBANK_SERVER}{massbank_id}'
    response = requests.get(request_url)
    spectrum_dict = response.json()
    mz, intensity = [], []
    for peak in spectrum_dict['spectrum'].split():
        peak_mz, peak_intensity = peak.split(':')
        mz.append(float(peak_mz))
        intensity.append(float(peak_intensity))
    precursor_mz = 0
    for metadata in spectrum_dict['metaData']:
        if metadata['name'] == 'precursor m/z':
            precursor_mz = float(metadata['value'])
            break
    source_link = (f'https://massbank.eu/MassBank/'
                   f'RecordDisplay.jsp?id={massbank_id}')
    return sus.MsmsSpectrum(usi, precursor_mz, 1, mz, intensity), source_link
