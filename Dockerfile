FROM mmtobservatory/mmtwfs:latest

MAINTAINER T. E. Pickering "te.pickering@gmail.com"

COPY . .

RUN pip install redis

RUN python setup.py develop

EXPOSE 8080

ENV WFSROOT /wfsdat

CMD ["wfssrv"]
