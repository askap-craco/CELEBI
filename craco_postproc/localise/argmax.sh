# Print the value in the second list at the index of the maximum of the
# first list.
# Adapted from https://stackoverflow.com/questions/62194066/optimally-finding-the-index-of-the-maximum-element-in-bash-array
# Usage: ./argmax.sh "1 2 3" "a b c"
awk 'BEGIN {
split(ARGV[1], a1);
split(ARGV[2], a2);
max=a1[1];
indx=1;
for (i in a1) {
	if (a1[i] > max) {
		indx = i;
		max = a1[i];
	}
}
print a2[indx]
}' "$@[1]" "$@[2]"
