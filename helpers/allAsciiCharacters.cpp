#include <iostream>
#include <cstdio>

int main(void) {
    
for(unsigned int c = '!'; c < ('!'+ 94); c++) {
    printf("%c", c);
}
    printf("\n");
        printf("\x1b[1;37;40m");
	for(unsigned int c = '!'; c < ('!'+ 94); c++) {
		printf("%c", c);
	}

    
        printf("\x1b[0m");
    printf("\n");
    
        printf("\x1b[2;37;40m");
    
for(unsigned int c = '!'; c < ('!'+ 94); c++) {
    printf("%c", c);
}
    printf("\x1b[0m");
    printf("\n");
	return 0;
}
