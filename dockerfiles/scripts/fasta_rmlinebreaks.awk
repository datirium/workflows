#!/usr/bin/gawk -f
{
	if(NR == 1) {
		printf("%s\n", $0);
	} else {
		if(substr($0,1,1) == ">") {
			printf("\n%s\n", $0);
		} else {
			printf("%s", $0);
		}
	}
}
END {
	printf("\n");
}
