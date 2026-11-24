FROM python:3.12-slim

RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH \
    PYTHONUNBUFFERED=1

WORKDIR /app

COPY --chown=user pyproject.toml README.md ./
COPY --chown=user src ./src
RUN pip install --no-cache-dir .

# Render/HF inyectan PORT; por defecto 7860 (puerto de Hugging Face Spaces)
CMD nanodelfy web --host 0.0.0.0 --port ${PORT:-7860} --no-browser
