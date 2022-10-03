To compile the all source code, you should first install the Boost C++ Libraries to your C++ compiler.

For Windows, please refer to the guide here: https://www.boost.org/doc/libs/1_77_0/more/getting_started/windows.html

For Linux and macOS, please refer to the guide here: https://www.boost.org/doc/libs/1_77_0/more/getting_started/unix-variants.html

On macOS devices, you can quickly install the Boost C++ Libraries using brew.

To install brew, run the below line in terminal
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

To install the Boost libraries with brew, run the below line terminal
brew install boost

For macOS and Linux:
After having Boost install, you can compile for getResults and makeIndex by redirecting to the folder of the source code in terminal. Then, you run "make". The executable files are compiled and placed in
the same folder.

For Windows:
Run the below lines to compile makeIndex.cpp and getResults.cpp respectively
g++ -std=c++17 makeIndex.cpp GridIndex.cpp -o makeIndex
g++ -std=c++17 getResults.cpp GridIndex.cpp -o getResults

Thank you, I wish and pray that you compile everything successfully.
If you have any problems compiling the source code, please tell me immediately.
