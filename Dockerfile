# SBONEST command-line image: docker run --rm -v "$PWD:/work" IMAGE check fit.json
FROM python:3.12-slim
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg \
    PIP_NO_CACHE_DIR=1 SBONEST_JOBS_DIR=/work/SB_JOBS
COPY . /src
RUN pip install -c /src/constraints-sideband.txt /src && rm -rf /src
WORKDIR /work
ENTRYPOINT ["sbonest"]
CMD ["--help"]
