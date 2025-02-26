FROM python:3.13

LABEL maintainer="te.pickering@gmail.com"

COPY . .

RUN apt-get update
RUN apt-get install -y git

RUN python -m pip install --upgrade pip setuptools setuptools_scm
RUN python -m pip install git+https://github.com/MMTObservatory/camsrv#egg=camsrv
RUN python -m pip install git+https://github.com/MMTObservatory/mmtwfs#egg=mmtwfs
RUN python -m pip install .

EXPOSE 8080

ENV WFSROOT /wfsdat

ENTRYPOINT ["wfssrv"]
