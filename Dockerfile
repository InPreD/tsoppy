# Stage 1: Builder
FROM python:3.14-slim-trixie AS builder

# Install git
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        git=1:2.47.3-0+deb13u1 \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Change to app dir and add mark as safe for git
WORKDIR /usr/src/app
RUN git config --global --add safe.directory /usr/src/app

# Copy pyproject.toml and source code
COPY .git/ ./.git/
COPY pyproject.toml ./
COPY src/ ./src/

# Install build dependencies and build the package as a wheel
RUN pip install --no-cache-dir --upgrade pip==26.2.1 \
    && pip install --no-cache-dir build==1.6.1 \
    && python -m build --wheel --outdir dist/ .

# Stage 2: Runtime
FROM python:3.14-slim-trixie AS runtime

WORKDIR /usr/src/app

# Copy only the built wheel and any other necessary runtime files
COPY --from=builder /usr/src/app/dist/*.whl .

# Install the built package from the wheel file
RUN pip install --no-cache-dir ./*.whl