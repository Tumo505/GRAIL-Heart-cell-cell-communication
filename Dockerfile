FROM python:3.10-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY requirements-dev.txt .
RUN pip install --no-cache-dir -r requirements-dev.txt

# Copy source code
COPY src/ ./src/
COPY setup.py .
COPY config.yaml .

# Install HeartMAP package
RUN pip install -e .[all]

# Expose port for API
EXPOSE 8000

# Default command
CMD ["python", "-m", "heartmap.api"]
