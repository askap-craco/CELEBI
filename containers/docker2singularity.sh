docker run -v /var/run/docker.sock:/var/run/docker.sock \
-v ./out:/output \
--privileged -t --rm \
quay.io/singularity/docker2singularity \
--name cracofunew-$(date +%y-%m-%d) cracofunew:latest