import sys, subprocess

# Trimming tools
def check_fastp_is_installed():
    try:
        subprocess.run(['fastp', '--version'], check=True)
        return True
    except subprocess.CalledProcessError:
        print("ERROR: Fastp is not installed or not in your PATH.")
        sys.exit(1)


def check_trimmomatic_is_installed():
    try:
        subprocess.run(['trimmomatic', '--version'], check=True)
        return True
    except subprocess.CalledProcessError:
        print("ERROR: Trimmomatic is not installed or not in your PATH.")
        sys.exit(1)


# Mapping tools
def check_STAR_is_installed():
    try:
        subprocess.run(['STAR', '--version'], check=True)
        return True
    except subprocess.CalledProcessError:
        print("ERROR: STAR is not installed or not in your PATH.")
        sys.exit(1)


# Assembly tools
def check_megahit_is_installed():
    try:
        subprocess.run(['megahit', '--version'], check=True)
        return True
    except subprocess.CalledProcessError:
        print("ERROR: Megahit is not installed or not in your PATH.")
        sys.exit(1)


def check_spades_is_installed():
    try:
        subprocess.run(['spades.py', '--version'], check=True)
        return True
    except subprocess.CalledProcessError:
        print("ERROR: SPAdes is not installed or not in your PATH.")
        sys.exit(1)


def check_rnaspades_is_installed():
    try:
        subprocess.run(['rnaspades.py', '--version'], check=True)
        return True
    except subprocess.CalledProcessError:
        print("ERROR: rnaspades is not installed or not in your PATH.")
        sys.exit(1)


def check_trinity_is_installed():
    try:
        subprocess.run(['Trinity', '--version'], check=True)
        return True
    except subprocess.CalledProcessError:
        print("ERROR: Trinity is not installed or not in your PATH.")
        sys.exit(1)

