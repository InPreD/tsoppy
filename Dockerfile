# Stage 1: Builder
FROM python:3.14-slim-trixie AS builder

WORKDIR /usr/src/app

# Copy pyproject.toml and source code
COPY pyproject.toml ./
COPY src/ ./src/

# Install build dependencies and build the package as a wheel
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir build \
    && python -m build --wheel --outdir dist/ .

# Stage 2: Runtime
FROM python:3.14-slim-trixie AS runtime

WORKDIR /usr/src/app

# Copy only the built wheel and any other necessary runtime files
COPY --from=builder /usr/src/app/dist/*.whl .

# Install the built package from the wheel file
RUN pip install --no-cache-dir ./*.whl