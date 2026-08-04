FROM python:3.11-slim

WORKDIR /app
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN python -m pip install --no-cache-dir --upgrade pip && \
    python -m pip install --no-cache-dir .

ENTRYPOINT ["reach"]
CMD ["--help"]
