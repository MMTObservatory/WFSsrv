FROM mmtobservatory/mmtwfs:latest

MAINTAINER T. E. Pickering "te.pickering@gmail.com"

COPY . .

RUN python -m pip install --upgrade pip
RUN python -m pip install git+https://github.com/MMTObservatory/camsrv.git#egg=camsrv
RUN python -m pip install git+https://github.com/MMTObservatory/cwfs.git#egg=cwfs
RUN python -m pip install -e .[all]

EXPOSE 8080

ENV WFSROOT /wfsdat

ENTRYPOINT ["wfssrv"]
