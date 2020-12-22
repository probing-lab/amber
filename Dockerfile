FROM python:3.8-slim-buster

WORKDIR /var
RUN apt-get update
RUN apt-get install -y git
RUN git clone https://github.com/probing-lab/amber.git 
WORKDIR /var/amber
RUN chmod +x /var/amber/amber
RUN pip install -r /var/amber/requirements.txt

ENTRYPOINT bash
