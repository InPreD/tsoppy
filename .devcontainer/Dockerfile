# Make hadolint available
FROM hadolint/hadolint:latest-debian AS hadolint

# Build on python devcontainer image
FROM mcr.microsoft.com/devcontainers/python:3.14-trixie

# Define build arguments for user configuration
ARG USERNAME=inpred
ARG USER_UID=1000
ARG USER_GID=$USER_UID
ARG FULLNAME="Your full name"
ARG MAIL="some@mail.address"

# Install packages and set locale
COPY --from=hadolint /bin/hadolint /bin/hadolint
RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
        silversearcher-ag=2.2.0+git20200805-1.2 \
        fish=4.0.2-1 \
    && echo "LC_ALL=en_US.UTF-8" > /etc/default/locale \
    && locale-gen en_US.UTF-8 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create the user
RUN groupadd --gid 1000 $USERNAME 2>/dev/null || true \
    && useradd --uid $USER_UID --gid 1000 -m $USERNAME -s /usr/bin/fish \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME
ENV PATH=$PATH:/home/$USERNAME/.local/bin
USER $USERNAME

# Configure git
RUN git config --global pull.rebase false && \
    git config --global user.email "${MAIL}" && \
    git config --global user.name "${FULLNAME}"
