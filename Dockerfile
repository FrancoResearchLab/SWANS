FROM pond:1.1

LABEL maintainer="Miguel Brown (brownm28@chop.edu)"
LABEL description="SWANS, best built with platform linux/amd64"

# Copy the local git repo into the image
COPY . /SWANS

# Make scripts in src/scripts executable
RUN chmod +x /SWANS/src/scripts/*

# Add scripts directory to PATH
ENV PATH="/SWANS/src/scripts:${PATH}"