## Windows: cl /LD /Ox /Qpar FPAMFilter.cpp
## Windows: "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Tools\MSVC\14.20.27508\bin\Hostx86\x64\cl.exe" /LD /Ox /Qpar FPAMFilter.cpp
all:
	g++ -shared -fPIC -O3 -o libFPAMFilter.so FPAMFilter.cpp