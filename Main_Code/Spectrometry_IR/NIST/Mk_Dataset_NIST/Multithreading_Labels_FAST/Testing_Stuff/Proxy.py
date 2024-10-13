import requests
import socks
import socket
import ssl

def get_smiles_from_name(compound_name):
    # Set up the SOCKS5 proxy
    socks.set_default_proxy(socks.SOCKS5, "199.102.104.70", 4145)
    socket.socket = socks.socksocket

    # Create a custom SSL context that doesn't verify certificates
    ssl_context = ssl.create_default_context()
    ssl_context.check_hostname = False
    ssl_context.verify_mode = ssl.CERT_NONE

    # Try the original name first
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/IsomericSMILES/JSON"
    response = requests.get(url, verify=False)
    if response.status_code == 200:
        data = response.json()
        return data['PropertyTable']['Properties'][0]['IsomericSMILES']

    # If not found, try removing parentheses and their contents
    modified_name = compound_name.split('(')[0].strip()
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{modified_name}/property/IsomericSMILES/JSON"
    response = requests.get(url, verify=False)
    if response.status_code == 200:
        data = response.json()
        return data['PropertyTable']['Properties'][0]['IsomericSMILES']

    # If still not found, try replacing spaces with hyphens
    modified_name = compound_name.replace(' ', '-')
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{modified_name}/property/IsomericSMILES/JSON"
    response = requests.get(url, verify=False)
    if response.status_code == 200:
        data = response.json()
        return data['PropertyTable']['Properties'][0]['IsomericSMILES']

    return None

print(get_smiles_from_name("Sulfur Dioxide"))
