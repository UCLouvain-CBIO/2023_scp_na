# Code for the manuscript: *Revisiting the thorny issue of missing values in single-cell proteomics*

The repository contains all the material required to reproduce the 
figures in the article:

> Vanderaa, Christophe, and Laurent Gatto. 2023. “Revisiting the 
Thorny Issue of Missing Values in Single-Cell Proteomics.” arXiv 
[q-bio.QM]. arXiv. http://arxiv.org/abs/2304.06654.

![](graphical_abstract.png)

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This image is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

## Content of the repository

- `make_figure[123].R`: 3 scripts that produce the article's figures, 
  each figure being generated by a separate script.
- `utils.R`: an R script with functions used by the scripts above.
- `install_dependencies.R`: an R script that installs all the required
  dependency packages. 
- `figs/`: folder containing the generated figures. 
- `Dockerfile`: the instructions used to build the Docker image. 

## Replicate the figures locally

Each figure is generated in its own script. We garantee long term 
reproducibility by relying on `Docker`. You must first install 
Docker on your local machine. Then, pull the image from the
[DockerHub
repository](https://hub.docker.com/repository/docker/cvanderaa/2023_scp_na_docker).

```
docker pull cvanderaa/2023_scp_na_docker
```

Then clone this repository on your local machine:

```
git clone git@github.com:UCLouvain-CBIO/2023_scp_na
cd 2023_scp_na/
```

Start an Rstudio session within a Docker container using:

```
docker run -e PASSWORD=bioc -p 8787:8787 -v `pwd`:/home/rstudio/ cvanderaa/2023_scp_na_docker:latest
```

Note you should use `%CD%` instead of `pwd` when using Windows. 

Open your browser and go to http://localhost:8787. The USER is
`rstudio` and the password is `bioc`. See the [DockerHub
repository](https://hub.docker.com/repository/docker/cvanderaa/2023_scp_na_docker)
for more detailed information on getting started with `Docker`.

## Note to future self

### Build the docker image

(I don't think I'll ever have to rebuild the image, but who knows). 

Build the image (make sure to `cd` in the `2023_scp_na` local repo):

```
docker build -t cvanderaa/2023_scp_na_docker .
```

When complete, push the new image to DockerHub:

```
docker push cvanderaa/2023_scp_na_docker
```

## Citation

To reuse the code for publications, please cite:

> Vanderaa, Christophe, and Laurent Gatto. 2023. “Revisiting the 
Thorny Issue of Missing Values in Single-Cell Proteomics.” arXiv 
[q-bio.QM]. arXiv. http://arxiv.org/abs/2304.06654.

## Licence

The code is provided under a permissive
[Artistic 2.0 license](https://opensource.org/license/artistic-license-2-0-php/).
