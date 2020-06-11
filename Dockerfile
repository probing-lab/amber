FROM python:3.7-slim

ADD . /var/amber
RUN chmod +x /var/amber/amber
RUN pip install -r /var/amber/requirements.txt
WORKDIR /var/amber