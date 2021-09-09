FROM archlinux:base-devel

# Update Arch
RUN pacman -Syu --noconfirm

# Install Git, sudo, python3, and pip
RUN pacman -S git sudo python python-pip --noconfirm

# Securely install yay using Git and makepkg on a separate user
RUN useradd builduser -m &&\
		passwd -d builduser &&\
		printf 'builduser ALL=(ALL) ALL\n' | tee -a /etc/sudoers &&\
		sudo -u builduser bash -c 'cd ~ && git clone https://aur.archlinux.org/yay.git && cd yay && makepkg -si --noconfirm'

# Install samtools
RUN sudo -u builduser yay -S samtools --noconfirm

# Install nimble
ADD . /nimble
RUN ["pip", "install", "./nimble"]

# Download the latest aligner version
RUN ["python3", "-m", "nimble", "download"]

ENTRYPOINT ["python3", "-m", "nimble"]
