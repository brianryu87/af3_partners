import urllib.request

TIMEOUT = 30


def get(url, accept=None):
    """GET a URL, return decoded body text. Raises urllib.error.URLError on failure."""
    req = urllib.request.Request(url)
    if accept:
        req.add_header("Accept", accept)
    with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
        return resp.read().decode("utf-8")
