# Make hadolint available
FROM hadolint/hadolint:latest-debian AS hadolint

# Build on python devcontainer image
FROM mcr.microsoft.com/devcontainers/python:3.14-trixie

# Define build arguments for user configuration
ARG USERNAME
ARG USER_UID
ARG FULLNAME
ARG MAIL

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

# Make /data directory
RUN mkdir /data \
    && chmod 777 /data

# Create the user
RUN groupadd --gid $USER_UID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_UID -m $USERNAME -s /usr/bin/fish\
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME
ENV PATH=$PATH:/home/$USERNAME/.local/bin
USER $USERNAME

# Configure git
RUN git config --global pull.rebase false && \
    git config --global user.email "${MAIL}" && \
    git config --global user.name "${FULLNAME}"
