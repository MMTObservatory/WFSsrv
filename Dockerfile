FROM python:3.13

LABEL maintainer="te.pickering@gmail.com"

RUN python -m pip install --upgrade pip
RUN python -m pip install -e .

EXPOSE 8080

ENV WFSROOT /wfsdat

ENTRYPOINT ["wfssrv"]
