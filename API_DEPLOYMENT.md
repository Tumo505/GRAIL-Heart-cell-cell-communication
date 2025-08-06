# HeartMAP API Deployment Guide

## Quick Start

### Local Development
1. Install dependencies:
```bash
pip install -r requirements-api.txt
```

2. Start the API server:
```bash
uvicorn api_server:app --host 0.0.0.0 --port 8001 --reload
```

3. Access the API:
- Interactive docs: http://localhost:8001/docs
- API root: http://localhost:8001/
- Health check: http://localhost:8001/health

### Docker Deployment

1. Build the Docker image:
```bash
docker build -t heartmap-api .
```

2. Run the container:
```bash
docker run -p 8001:8000 heartmap-api
```

3. Access the API at http://localhost:8001

### Production Deployment

#### Using Docker Compose

Create `docker-compose.yml`:
```yaml
version: '3.8'
services:
  heartmap-api:
    build: .
    ports:
      - "8001:8000"
    environment:
      - HEARTMAP_DATA_DIR=/app/data
    volumes:
      - ./data:/app/data
      - ./uploads:/app/uploads
    restart: unless-stopped
```

Run with:
```bash
docker-compose up -d
```

#### Cloud Deployment Options

1. **AWS ECS/Fargate**
2. **Google Cloud Run**
3. **Azure Container Instances**
4. **Heroku**
5. **Railway**
6. **DigitalOcean App Platform**

## API Endpoints

### Core Endpoints
- `GET /` - API information
- `GET /health` - Health check
- `GET /models` - Available analysis models
- `POST /analyze` - Run analysis on uploaded data

### Configuration
- `GET /config` - Get current configuration
- `POST /config` - Update configuration

## Usage Examples

### Upload and Analyze Data

```python
import requests

# Upload file and analyze
with open('heart_data.h5ad', 'rb') as f:
    files = {'file': f}
    data = {
        'analysis_type': 'comprehensive',
        'output_format': 'json'
    }
    response = requests.post('http://localhost:8001/analyze', files=files, data=data)
    results = response.json()
```

### Using curl

```bash
curl -X POST "http://localhost:8001/analyze" \
  -F "file=@heart_data.h5ad" \
  -F "analysis_type=basic"
```

## Environment Variables

- `HEARTMAP_DATA_DIR` - Directory for data storage
- `HEARTMAP_CONFIG_PATH` - Path to configuration file
- `HEARTMAP_LOG_LEVEL` - Logging level (INFO, DEBUG, etc.)

## Security Considerations

- Enable HTTPS in production
- Implement authentication/authorization
- Set up rate limiting
- Configure CORS properly
- Use environment variables for secrets

## Monitoring

- Health check endpoint: `/health`
- Prometheus metrics (can be added)
- Log aggregation
- Application performance monitoring

## Scaling

- Use multiple workers with gunicorn
- Implement horizontal scaling with load balancers
- Consider async processing for large files
- Use Redis/Celery for background tasks
