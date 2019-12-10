FROM mmtobservatory/mmtwfs:latest

MAINTAINER T. E. Pickering "te.pickering@gmail.com"

COPY . .

RUN pip install -r requirements.txt
RUN python setup.py install

EXPOSE 8080

ENV WFSROOT /wfsdat

CMD ["wfssrv"]
