FROM ubuntu:22.04

# Prevent interactive prompts during apt installations
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies: python3, git, build essentials, and zlib
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    git \
    make \
    gcc \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone and build TRF-mod from source
RUN git clone https://github.com/lh3/TRF-mod.git /tmp/TRF-mod && \
    cd /tmp/TRF-mod && \
    make -f compile.mak && \
    cp trf-mod /usr/local/bin/ && \
    rm -rf /tmp/TRF-mod

# Set up the application directory
WORKDIR /app

# Copy pipeline scripts into the container
COPY *.py /app/
COPY run_wrapper.sh /app/

# Ensure the bash script is executable
RUN chmod +x /app/run_wrapper.sh

# Set the entrypoint to the wrapper script
ENTRYPOINT ["/app/run_wrapper.sh"]

# Default command if no arguments are provided (prints usage)
CMD ["--help"]
