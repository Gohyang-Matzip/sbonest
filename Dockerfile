# SBONEST command-line image: docker run --rm -v "$PWD:/work" IMAGE check fit.json
# optimalcontrol-nmr compiles a Rust extension from source, so a full image builds
# the environment and the slim runtime image only receives the result.
FROM python:3.12 AS build
ENV PIP_NO_CACHE_DIR=1
COPY . /src
RUN python -m venv /opt/sbonest \
    && /opt/sbonest/bin/pip install -c /src/constraints-sideband.txt /src

FROM python:3.12-slim
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg \
    SBONEST_JOBS_DIR=/work/SB_JOBS PATH=/opt/sbonest/bin:$PATH
COPY --from=build /opt/sbonest /opt/sbonest
WORKDIR /work
ENTRYPOINT ["sbonest"]
CMD ["--help"]
