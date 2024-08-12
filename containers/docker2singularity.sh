docker run -v /var/run/docker.sock:/var/run/docker.sock \
-v .:/output \
--privileged -t --rm \
quay.io/singularity/docker2singularity \
-name cracofunew cracofunew:latest