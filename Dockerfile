FROM mmtobservatory/mmtwfs:latest

MAINTAINER T. E. Pickering "te.pickering@gmail.com"

RUN pip install redis

COPY . .

RUN python setup.py develop

EXPOSE 8080

ENV WFSROOT /wfsdat

ENTRYPOINT ["wfssrv"]
