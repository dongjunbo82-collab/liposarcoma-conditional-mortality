"""Download only public reference inputs; never send SEER records over the network."""
from pathlib import Path
import argparse,os,time,urllib.request,urllib.parse
from lxml import etree as E
PARAMETERS = {'D76': {'B_1': ['D76.V1-level1'], 'B_2': ['D76.V7'], 'B_3': ['D76.V51'], 'B_4': ['*None*'], 'B_5': ['*None*'], 'F_D76.V1': ['2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017'], 'F_D76.V10': ['*All*'], 'F_D76.V2': ['I00', 'I01', 'I02', 'I05', 'I06', 'I07', 'I08', 'I09', 'I10', 'I11', 'I12', 'I13', 'I20', 'I21', 'I22', 'I24', 'I25', 'I26', 'I27', 'I28', 'I30', 'I31', 'I33', 'I34', 'I35', 'I36', 'I37', 'I38', 'I40', 'I42', 'I44', 'I45', 'I46', 'I47', 'I48', 'I49', 'I50', 'I51', 'I60', 'I61', 'I62', 'I63', 'I64', 'I67', 'I69', 'I70', 'I71', 'I72', 'I73', 'I74', 'I77', 'I78'], 'F_D76.V25': ['*All*'], 'F_D76.V27': ['*All*'], 'F_D76.V9': ['*All*'], 'M_1': ['D76.M1'], 'M_2': ['D76.M2'], 'M_3': ['D76.M3'], 'O_V10_fmode': ['freg'], 'O_V1_fmode': ['freg'], 'O_V25_fmode': ['freg'], 'O_V27_fmode': ['freg'], 'O_V2_fmode': ['freg'], 'O_V9_fmode': ['freg'], 'O_age': ['D76.V51'], 'O_javascript': ['on'], 'O_location': ['D76.V9'], 'O_oc-sect1-request': ['close'], 'O_precision': ['5'], 'O_rate_per': ['100000'], 'O_show_suppressed': ['true'], 'O_show_totals': ['false'], 'O_show_zeros': ['true'], 'O_timeout': ['600'], 'O_ucd': ['D76.V2'], 'O_urban': ['D76.V19'], 'VM_D76.M6_D76.V10': [''], 'VM_D76.M6_D76.V17': ['*All*'], 'VM_D76.M6_D76.V1_S': ['*All*'], 'VM_D76.M6_D76.V7': ['*All*'], 'VM_D76.M6_D76.V8': ['*All*'], 'V_D76.V1': [''], 'V_D76.V10': [''], 'V_D76.V11': ['*All*'], 'V_D76.V12': ['*All*'], 'V_D76.V17': ['*All*'], 'V_D76.V19': ['*All*'], 'V_D76.V2': ['*All*'], 'V_D76.V20': ['*All*'], 'V_D76.V21': ['*All*'], 'V_D76.V22': ['*All*'], 'V_D76.V23': ['*All*'], 'V_D76.V24': ['*All*'], 'V_D76.V25': [''], 'V_D76.V27': [''], 'V_D76.V4': ['*All*'], 'V_D76.V5': ['*All*'], 'V_D76.V51': ['*All*'], 'V_D76.V52': ['*All*'], 'V_D76.V6': ['00'], 'V_D76.V7': ['*All*'], 'V_D76.V8': ['*All*'], 'V_D76.V9': [''], 'action-Send': ['Send'], 'dataset_code': ['D76'], 'dataset_label': ['Underlying Cause of Death, 1999-2016'], 'dataset_vintage': ['2016'], 'stage': ['request'], 'accept_datause_restrictions': ['true'], 'O_aar_pop': ['0000'], 'O_aar_enable': [''], 'O_aar': ['aar_none'], 'O_aar_CI': [''], 'O_aar_SE': ['']}, 'D158': {'query': [''], 'affiliate': ['cdc-main'], 'dataset_code': ['D158'], 'dataset_label': ['Underlying Cause of Death, 2018-2024, Single Race'], 'dataset_vintage': ['2024'], 'stage': ['request'], 'O_javascript': ['off'], 'M_1': ['D158.M1'], 'M_2': ['D158.M2'], 'M_3': ['D158.M3'], 'B_1': ['D158.V1-level1'], 'B_2': ['D158.V7'], 'B_3': ['D158.V51'], 'B_4': ['*None*'], 'B_5': ['*None*'], 'O_oc-sect1-request': [''], 'O_rate_per': ['100000'], 'VM_D158.M6_D158.V1_S': ['*All*'], 'VM_D158.M6_D158.V7': ['*All*'], 'VM_D158.M6_D158.V17': ['*All*'], 'VM_D158.M6_D158.V42': ['*All*'], 'VM_D158.M6_D158.V10': [''], 'O_location': ['D158.V9'], 'O_V9_fmode': ['freg'], 'V_D158.V9': [''], 'F_D158.V9': ['*All*'], 'O_V10_fmode': ['freg'], 'V_D158.V10': [''], 'F_D158.V10': ['*All*'], 'O_V27_fmode': ['freg'], 'V_D158.V27': [''], 'F_D158.V27': ['*All*'], 'O_V30_fmode': ['freg'], 'V_D158.V30': [''], 'F_D158.V30': ['*All*'], 'O_V31_fmode': ['freg'], 'V_D158.V31': [''], 'F_D158.V31': ['*All*'], 'O_urban': ['D158.V19'], 'V_D158.V19': ['*All*'], 'V_D158.V11': ['*All*'], 'V_D158.V18': ['*All*'], 'O_age': ['D158.V51'], 'V_D158.V5': ['*All*'], 'V_D158.V51': ['*All*'], 'V_D158.V52': ['*All*'], 'V_D158.V6': ['00'], 'V_D158.V7': ['*All*'], 'V_D158.V17': ['*All*'], 'O_race': ['D158.V42'], 'V_D158.V42': ['*All*'], 'V_D158.V43': ['*All*'], 'V_D158.V44': ['*All*'], 'V_D158.V45': ['*All*'], 'O_V1_fmode': ['freg'], 'V_D158.V1': [''], 'F_D158.V1': ['2018', '2019', '2020', '2021', '2022', '2023'], 'V_D158.V24': ['*All*'], 'V_D158.V20': ['*All*'], 'V_D158.V21': ['*All*'], 'O_ucd': ['D158.V2'], 'O_V2_fmode': ['freg'], 'V_D158.V2': ['*All*'], 'F_D158.V2': ['I00', 'I01', 'I02', 'I05', 'I06', 'I07', 'I08', 'I09', 'I10', 'I11', 'I12', 'I13', 'I20', 'I21', 'I22', 'I24', 'I25', 'I26', 'I27', 'I28', 'I30', 'I31', 'I33', 'I34', 'I35', 'I36', 'I37', 'I38', 'I40', 'I42', 'I44', 'I45', 'I46', 'I47', 'I48', 'I49', 'I50', 'I51', 'I60', 'I61', 'I62', 'I63', 'I64', 'I67', 'I69', 'I70', 'I71', 'I72', 'I73', 'I74', 'I77', 'I78'], 'V_D158.V4': ['*All*'], 'V_D158.V12': ['*All*'], 'V_D158.V22': ['*All*'], 'V_D158.V23': ['*All*'], 'O_V25_fmode': ['freg'], 'V_D158.V25': [''], 'F_D158.V25': ['*All*'], 'O_export-format': ['xls'], 'O_show_totals': ['false'], 'O_precision': ['5'], 'O_timeout': ['600'], 'O_show_zeros': ['true'], 'accept_datause_restrictions': ['true'], 'O_aar_pop': ['0000'], 'O_aar_enable': [''], 'O_aar': ['aar_none'], 'O_aar_CI': [''], 'O_aar_SE': ['']}}

def build_request(database, ten_year=False):
    values={k:list(v) for k,v in PARAMETERS[database].items()}
    if ten_year:
        values['B_3']=[database+'.V5'];values['O_age']=[database+'.V5']
    root=E.Element('request-parameters')
    for k,vs in values.items():
        p=E.SubElement(root,'parameter');E.SubElement(p,'name').text=k
        for v in vs:E.SubElement(p,'value').text=v
    return E.tostring(root)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--seer-only',action='store_true');args=parser.parse_args()
    base=Path(os.environ.get('LPS_OUTPUT_DIR','outputs/population_study'))
    target=base/'survivorship';target.mkdir(parents=True,exist_ok=True)
    request=urllib.request.Request('https://seer.cancer.gov/codrecode/1969_d03012018/index.txt',headers={'User-Agent':'Research-reproducibility-script/1.0'})
    with urllib.request.urlopen(request,timeout=60) as response:content=response.read()
    (target/'seer_cod_recode.txt').write_bytes(content)
    if args.seer_only:return
    target=base/'revision/wonder_20261016';target.mkdir(parents=True,exist_ok=True)
    for i,(db,span,ten) in enumerate([('D76','2009_2017',False),('D76','2009_2017',True),('D158','2018_2023',False),('D158','2018_2023',True)]):
        if i:time.sleep(16)
        xml=build_request(db,ten)
        stem=db+'_'+span+('_ten_year' if ten else '')
        (target/(stem+'_request.xml')).write_bytes(xml)
        payload=urllib.parse.urlencode({'request_xml':xml.decode(),'accept_datause_restrictions':'true'}).encode()
        req=urllib.request.Request('https://wonder.cdc.gov/controller/datarequest/'+db,data=payload)
        with urllib.request.urlopen(req,timeout=60) as response:content=response.read()
        root=E.fromstring(content)
        if not root.findall('.//data-table/r'):raise RuntimeError('No reference rows returned; review the current CDC service and its access conditions.')
        (target/(stem+'_response.xml')).write_bytes(content)
    print('Public reference responses downloaded; next run analysis/validate_wonder_reference.py.')

if __name__=='__main__':main()
