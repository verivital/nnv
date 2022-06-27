FROM ubuntu:20.04

RUN apt-get update && \
    apt-get install -y python3.8 python3.8-dev python3-pip && \
    apt-get clean

RUN apt-get install -qy curl unzip

RUN curl -fL https://github.com/mldiego/GoTube/archive/refs/heads/main.zip > gotube.zip

RUN unzip gotube.zip

WORKDIR GoTube-main

RUN chmod +x paper_submission_compare.sh

RUN pip3 install -r requirements.txt

RUN chmod +x paper_submission_short.sh

# Clear previous results from GitHub
RUN mv /GoTube-main/logged /GoTube-main/orig_logged
RUN mv /GoTube-main/saved_outputs /GoTube-main/orig_saved_outputs

# Long version, about 2-3 hours
# RUN ./paper_submission_compare.sh

# Short version (max 1 hour)
# RUN ./paper_submission_short.sh
