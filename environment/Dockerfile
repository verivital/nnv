FROM registry.codeocean.com/codeocean/matlab:2017b-ubuntu16.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

COPY postInstall /
RUN /postInstall
