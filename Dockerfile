FROM python:3.10-slim

RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/IGMT

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

# PLINK
RUN mkdir -p /opt/software/plink1.9 && \
    cd /opt/software/plink1.9 && \
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip && \
    unzip plink_linux_x86_64_20250819.zip && \
    rm plink_linux_x86_64_20250819.zip && \
    chmod +x plink

# EMMAX
RUN mkdir -p /opt/software/EMMAX && \
    cd /opt/software/EMMAX && \
    wget https://csg.sph.umich.edu/kang/emmax/download/emmax-intel-binary-20120210.tar.gz && \
    tar -xzf emmax-intel-binary-20120210.tar.gz && \
    rm emmax-intel-binary-20120210.tar.gz && \
    chmod +x emmax-intel64 emmax-kin-intel64

COPY . /opt/IGMT

RUN chmod +x /opt/IGMT/IGMT.py

ENV PATH="/opt/IGMT:${PATH}"
ENV PYTHONPATH="/opt/IGMT/src:${PYTHONPATH}"

ENTRYPOINT ["python3", "/opt/IGMT/IGMT.py"]
CMD ["-h"]
