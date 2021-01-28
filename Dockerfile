FROM python:3.8-slim-buster

WORKDIR /var
COPY . /var/amber
WORKDIR /var/amber
RUN chmod +x /var/amber/amber
RUN pip install -r /var/amber/requirements.txt

ENTRYPOINT bash
