import re

def ftp_to_https(ftp_url):
    """Convert a URL from FTP to HTTPS protocol."""
    https_url = re.sub(r'^ftp:\/\/(.+)$', 'https://\\1', ftp_url)
    print('FTP {}\n-> HTTPS {}'.format(ftp_url, https_url))
    return https_url
