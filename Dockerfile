FROM mmtobservatory/mmtwfs:latest

MAINTAINER T. E. Pickering "te.pickering@gmail.com"

COPY . .

RUN pip install git+https://github.com/MMTObservatory/camsrv.git#egg=camsrv
RUN pip install -e [all]

EXPOSE 8080

ENV WFSROOT /wfsdat

ENTRYPOINT ["wfssrv"]
