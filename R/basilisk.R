pythonEnv <- BasiliskEnvironment(
    "env1",
    pkgname = "oob",
    channels = c("conda-forge", "bioconda"),
    pip = c("trimap==1.0.15",
        "igraph==0.11.6",
                "leidenalg==0.10.2"),
    packages = c(
        "numpy==1.21.0"
    )
)
