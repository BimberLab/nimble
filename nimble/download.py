import requests
from nimble.utils import get_exec_name_from_platform

# Given the name of a release, download the platform-specific executable from that release.
# Given no name, default to the most recent release.
def download(release):
    exec_name = get_exec_name_from_platform()
    print("Downloading " + exec_name)

    # Get a specific release, if one was provided. Otherwise get the most recent
    url = ""
    if len(release) == 1:
        url = (
            "https://github.com/BimberLab/nimble-aligner/releases/download/"
            + release[0]
            + "/"
            + exec_name
        )
    else:
        url = (
            "https://github.com/BimberLab/nimble-aligner/releases/latest/download/"
            + exec_name
        )

    # Download aligner
    r = requests.get(url)

    # The executable will be placed in the python package directory
    aligner_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "aligner")
    print("Aligner download path: " + aligner_path)

    # If we successfully downloaded it, write the file to disk and give it execution permissions if necessary
    if r.status_code == 200:
        with open(aligner_path, "wb") as f:
            f.write(r.content)
        if sys.platform == "linux" or sys.platform == "darwin": # If we're on Unix, specify permission to write
            st = os.stat(aligner_path)
            os.chmod(aligner_path, st.st_mode | stat.S_IEXEC | stat.S_IXOTH)
    else:
        print("Error -- could not download aligner, status code " + str(r.status_code))
        sys.exit()