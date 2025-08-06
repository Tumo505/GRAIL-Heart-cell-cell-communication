FROM python:3.10-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY requirements-api.txt .
RUN pip install --no-cache-dir -r requirements-api.txt

# Copy the API server
COPY api_server.py .

# Create a directory for data uploads
RUN mkdir -p /app/uploads

# Expose port
EXPOSE 8000

# Set environment variables
ENV PYTHONPATH=/app
ENV HEARTMAP_DATA_DIR=/app/data

# Run the application
CMD ["python", "api_server.py"]

# Expose port for API
EXPOSE 8000

# Default command
CMD ["python", "-m", "heartmap.api"]
