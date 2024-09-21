pythonEnv <- BasiliskEnvironment(
    "env1",
    pkgname = "oob",
    channels = c("conda-forge", "bioconda"),
    pip = c("igraph==0.11.6",
        "leidenalg==0.10.2",
        "numpy==1.26.4"),
    packages = c(
        "trimap==1.0.15"
    )
)
