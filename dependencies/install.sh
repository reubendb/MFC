#!/usr/bin/env bash

# Logo generated here - and then modified - https://patorjk.com/software/taag/#p=moreopts&c=echo&f=Isometric2&t=MFC.
function print_logo() {
    echo -en "\033[0;34m"
    echo "      ___           ___           ___     "
    echo "     /\\  \\         /\\__\\         /\\__\\    "
    echo "    |##\\  \\       /#/ _/_       /#/  /    "
    echo "    |#|#\\  \\     /#/ /\\__\\     /#/  /     "
    echo "  __|#|\#\\  \\   /#/ /#/  /    /#/  /  ___ "
    echo " /####|_\#\\__\\ /#/_/#/  /    /#/__/  /\\__\\"
    echo " \\#\\~~\\  \\/__/ \\#\\/#/  /     \\#\\  \\ /#/  /"
    echo "  \\#\\  \\        \\##/__/       \\#\\  /#/  / "
    echo "   \\#\\  \\        \\#\\  \\        \\#\\/#/  /  "
    echo "    \\#\\__\\        \\#\\__\\        \\##/  /   "
    echo "     \\/__/         \\/__/         \\/__/    "
    echo -e "\033[0m"
}

print_logo

echo -e "${COLORS[ORANGE]}\nPlease run ./install.sh <-h|--help> for more information.\n${COLORS[NONE]}"

# Helper Functions

function print_centered_text() { printf "%*s\n" $(((${#1}+$(tput cols))/2)) "$1" ; }
function print_line_of_char()  { printf "%$(tput cols)s\n" | tr " " "$1"         ; }
function print_bounded_line()  { printf "|| %-$(($(tput cols)-6))s ||" "$1"      ; }
function clear_line()          { echo -en "\r\033[2K"                            ; }

# Error Handling

function on_exit() {
    exit_code=$?

    if [[ exit_code -ne "0" ]]; then
        echo -e "\033[0;31m"
        print_line_of_char "="
        print_bounded_line "Fatal Error"
        print_bounded_line "---> Exit Code ($exit_code)."
        print_bounded_line "---> Take a look at the output above."
        print_bounded_line "---> Some stages output logs to log/."
        print_bounded_line "---> Please ensure your build environment is adequate (view README.md)."
        print_line_of_char "="
        echo -e "\033[0m"
        exit $exit_code
    fi
}

trap on_exit EXIT ; set -e ; set -o pipefail ; set -o errtrace

# Check Shell Compatibility

SEHLL_TEST_COMMANDS=("declare -A __shell_test_dict__")

i=1
for shell_test_command in "${SEHLL_TEST_COMMANDS[@]}"; do
    clear_line
    echo -en "\r|--> ($i/${#SEHLL_TEST_COMMANDS[@]}) Checking whether \"$shell_test_command\" runs..."
    if ! ($shell_test_command); then
        clear_line
        echo -e "\r|--> \033[0;31mFailed to run \"$shell_test_command\". Please check your shell is supported and up to date\033[0m."
        echo -e "\r|--> \033[0;33mNotably, macOS is known to run older versions of Bash\033[0m."
        exit 1
    fi
    i=$((i+1))
done

clear_line
echo -e "\r|--> \033[0;32mSuccessfully ran shell compatibility checks\033[0m."

# Start of install.sh

declare -A DOWNLOADED_DEPENDENCIES COLORS            DEFAULT_UTILITIES \
           CHECK_UTILITIES         DIRECTORIES       USED_UTILITIES

DOWNLOADED_DEPENDENCIES=( ["FFTW3"]="http://www.fftw.org/fftw-3.3.10.tar.gz" )

COLORS=( [RED]="\033[0;31m"    [GREEN]="\033[0;32m" \
         [ORANGE]="\033[0;33m" [NONE]="\033[0m"     )

DEFAULT_UTILITIES=( [GIT]="git"         [MAKE]="make"                     \
                    [WGET]="wget|curl"  [TAR]="tar"                       \
                    [PYTHON3]="python3" [PYTHON3-CONFIG]="python3-config" \
                    [CC]="mpicc"        [CPPC]="mpicxx"                   \
                    [FC77]="mpif77"     [FC90]="mpif90"                   )

DOTFILE_FILENAMES=( ".bashrc" ".bash_profile" ".bash_login" ".profile" \
                    ".zshrc"  ".cshrc"        ".zprofile"              )

DIRECTORIES=( [DEPENDENCIES]="$(pwd)" [SRC]="$(pwd)/src" \
              [BUILD]="$(pwd)/build"  [LOG]="$(pwd)/log" )

SPIN_ANIMATION=("-" "\\" "|" "/")

SLEEP_DURATION=0.1

JOB_COUNT=1

CFLAGS=""
CPPFLAGS=""

if [ ! -f ./install.sh ]; then
    echo -e "${COLORS[RED]}Error: You must call this script from within the 'dependencies' folder.\n${COLORS[NONE]}"
    exit 1
fi

function show_help() {
    echo -e "Usage: ./install.sh ${COLORS[ORANGE]}[options]${COLORS[NONE]}                                                                                  "
    echo "Options:                                                                                                                                          "
    echo -e "  ${COLORS[ORANGE]}-h,       --help                   ${COLORS[NONE]} Prints this message and exits.                                           "
    echo -e "  ${COLORS[ORANGE]}-c,       --clean                  ${COLORS[NONE]} Cleans.                                                                  "
    echo -e "  ${COLORS[ORANGE]}-j    [N] --jobs                [N]${COLORS[NONE]} Allows for ${COLORS[ORANGE]}'N'${COLORS[NONE]} concurrent jobs.          "
    echo -e "  ${COLORS[ORANGE]}-cc   [X] --c-compiler          [X]${COLORS[NONE]} Uses ${COLORS[ORANGE]}'X'${COLORS[NONE]} as the C compiler.              "
    echo -e "  ${COLORS[ORANGE]}-cppc [X] --cpp-compiler        [X]${COLORS[NONE]} Uses ${COLORS[ORANGE]}'X'${COLORS[NONE]} as the CPP compiler.            "
    echo -e "  ${COLORS[ORANGE]}-fc77 [X] --fortran-compiler-77 [X]${COLORS[NONE]} Uses ${COLORS[ORANGE]}'X'${COLORS[NONE]} as the Fortran 77 compiler.     "
    echo -e "  ${COLORS[ORANGE]}-fc90 [X] --fortran-compiler-90 [X]${COLORS[NONE]} Uses ${COLORS[ORANGE]}'X'${COLORS[NONE]} as the Fortran 90 compiler.     "
    echo -e "  ${COLORS[ORANGE]}-cf   [X] --c-flags             [X]${COLORS[NONE]} Uses ${COLORS[ORANGE]}'X'${COLORS[NONE]} as the C compiler's flags.      "
    echo -e "  ${COLORS[ORANGE]}-cppf [X] --cpp-flags           [X]${COLORS[NONE]} Uses ${COLORS[ORANGE]}'X'${COLORS[NONE]} as the CPP compiler's flags.    "
    echo -e "  ${COLORS[ORANGE]}-ff   [X] --fortran-flags       [X]${COLORS[NONE]} Uses ${COLORS[ORANGE]}'X'${COLORS[NONE]} as the Fortran compiler's flags."

    echo -e "\nProvided courtesy of MFC (https://github.com/MFlowCode/MFC)."
}

function show_command_running() {
    base_string="$1" ; shift

    # A visial cue if "$@" fails right away
    echo -en "$base_string ${SPIN_ANIMATION[0]}"

    SECONDS=0
    "$@" &

    pid=$!
    while ps -p $pid &>/dev/null; do
        for ((i=0;i<${#SPIN_ANIMATION[@]};i++)); do
            echo -en "\r$base_string ${SPIN_ANIMATION[$i]} ($(("$SECONDS" / 60))m $(("$SECONDS" % 60))s)"
            sleep 0.25
        done
    done

    # Once the command has finished running, we run `wait` to get its exit code
    wait $pid
}

function log_command() {
    log_filepath="$1" ; shift
    "$@" >> "$log_filepath" 2>&1
}

while [[ $# -gt 0 ]]; do
    option="$1"
    case $option in
        -h   |--help               ) show_help                  ; exit  0       ;;
        -cc  |--c-compiler         ) CHECK_UTILITIES[CC]="$2"   ; shift ; shift ;;
        -fc77|--fortran-compiler-77) CHECK_UTILITIES[FC77]="$2" ; shift ; shift ;;
        -fc90|--fortran-compiler-90) CHECK_UTILITIES[FC90]="$2" ; shift ; shift ;;
        -cppc|--cpp-compiler       ) CHECK_UTILITIES[CPPC]="$2" ; shift ; shift ;;
        -j   |--jobs               ) JOB_COUNT="$2"             ; shift ; shift ;;
        -cf  |--c-flags            ) CFLAGS="$2"                ; shift ; shift ;;
        -cppf|--cpp-flags          ) CPPFLAGS="$2"              ; shift ; shift ;;
        -c   |--clean              )
            rm -rf "${DIRECTORIES[LOG]}" "${DIRECTORIES[BUILD]}"
            for dependency_name in "${!DOWNLOADED_DEPENDENCIES[@]}"; do
                rm -rf "${DIRECTORIES[SRC]}:?/$dependency_name"
            done
            exit 0
            ;;
        *                          )
            echo -e "${COLORS[RED]}Error: Unknown command line option \"$option\".${COLORS[NONE]}"
            echo "Please run ./install.sh <-h|--help> for more information."
            exit 1
            ;; 
    esac
done

for utility_name in "${!DEFAULT_UTILITIES[@]}"; do
    if [[ ! -v "CHECK_UTILITIES[$utility_name]" ]] ; then
        CHECK_UTILITIES[$utility_name]=${DEFAULT_UTILITIES[$utility_name]}
    fi
done

# Create required directories
mkdir -p "${DIRECTORIES[SRC]}" "${DIRECTORIES[BUILD]}" \
         "${DIRECTORIES[LOG]}"

function check_command_existance() {
    i=1 ; for utility_name in "${!CHECK_UTILITIES[@]}"; do
        IFS="|" read -r -a utility_alternatives <<< "${CHECK_UTILITIES[$utility_name]}"
        
        for utility_alternative in "${utility_alternatives[@]}"; do
            if command -v "$utility_alternative" &> /dev/null; then
                USED_UTILITIES[$utility_name]="$utility_alternative"        
                break
            fi
        done

        line_color="${COLORS[GREEN]}" ; status_text=" FOUND "
        using_text="${USED_UTILITIES[$utility_name]}" 

        if [[ "${USED_UTILITIES[$utility_name]}" == "" ]]; then
            line_color="${COLORS[RED]}"
            status_text="MISSING"
            using_text="?"
        fi

        printf "\r%5s Looking for %-40s: $line_color%7s${COLORS[NONE]} (Using %s)"     \
            "$i/${#CHECK_UTILITIES[@]}"                                                \
            "$utility_name ($(echo "${utility_alternatives[*]}" | sed s/\ /\ or\ /g))" \
            "$status_text" "$using_text"
        
        if [[ "${USED_UTILITIES[$utility_name]}" == "" ]]; then
            echo -e "\n|--> ${COLORS[RED]}Failed to find one of the above utilities.\n${COLORS[NONE]}"
            exit 1
        fi

        i=$((i+1))

        sleep $SLEEP_DURATION
    done

    clear_line
    echo -e "\r|--> ${COLORS[GREEN]}Found all required/selected command line utilities${COLORS[NONE]}:" \
            "$(echo "${USED_UTILITIES[*]}" | sed s/\ /,\ /g | rev | sed s/\ ,/\ dna\ ,/ | rev)"
}

check_command_existance

# If running on the Expanse supercomputer
if [[ $(env | grep -i 'expanse' | wc -c) -ne 0 ]]; then
    module load numactl
fi

show_command_running "|--> (1/2) Fetching Submodules... "    \
                     git submodule update --init --recursive
clear_line

cd "${DIRECTORIES[SRC]}"
    i=1 ; for dependency_name in "${!DOWNLOADED_DEPENDENCIES[@]}"; do
        dependency_link=${DOWNLOADED_DEPENDENCIES[$dependency_name]}

        base_string="|--> (2/2) Fetching Archives - ($i/${#DOWNLOADED_DEPENDENCIES[@]}) $dependency_name"

        # If we haven't downloaded it before (the directory named $name doesn't exist)
        if [ ! -d "$(pwd)/$dependency_name" ]; then
            archive_filename="$dependency_name.tar.gz"
            
            case ${USED_UTILITIES[WGET]} in
                "wget") show_command_running "$base_string: Downloading Archive..."                  \
                                             wget -O "$archive_filename" -q "$dependency_link"
                    ;;
                "curl") show_command_running "$base_string: Downloading Archive..."                  \
                                             curl -s "$dependency_link" --output "$archive_filename"
                    ;;
            esac

            clear_line ; mkdir -p "$dependency_name"

            cd "$dependency_name"
                show_command_running "$base_string: Uncompressing Downloaded Archive..." \
                                     tar --strip-components 1 -xf "../$archive_filename"
            cd ..

            rm "$archive_filename"
        fi

        i=$((i+1))
    done
cd ..

clear_line

N_DEPENDENCIES=$(find "${DIRECTORIES[SRC]}" -mindepth 1 -maxdepth 1 -type d | wc -l)

echo -en "\r|--> ${COLORS[GREEN]}Successfully fetched the source code of all $N_DEPENDENCIES dependencies${COLORS[NONE]}: "

i=1 ; for entry in "${DIRECTORIES[SRC]}"/*; do
    if [ -d "$entry" ]; then
        if [ "$i" -eq "$((N_DEPENDENCIES))" ]; then
            echo -n ", and "
        elif [ "$i" -ne "1" ]; then
            echo -n ", "
        fi

        echo -n "$(basename "$entry")"

        if [ "$i" -eq "$N_DEPENDENCIES" ]; then
            echo -en ".\n"
        fi 
        
        i=$((i+1))
    fi
done


# 2) Add "build/<lib, bin, include, ...>" to the system's search path
#
# We append the export USED_UTILITIES to every dotfile we can find because it's
# tricky to know which ones will be executed and when, especially if the
# user has multiple shells installed at the same time. There are also
# some complexities with login and non-login shell sessions.
#
# The code we add to the dotfiles has include guards to prevent
# errors and redundancy. Also, we only append to the dotfiles if
# the header guard isn't set to prevent duplicates.

found_dotfile_count=0
found_dotfile_list_string=""
if [[ -z "${MFC_ENV_SH_HEADER_GUARD}" ]]; then 
    export_cmds_0="export MFC_ENV_SH_HEADER_GUARD=\"SET\""
    export_cmds_1="export LD_LIBRARY_PATH=\"\$LD_LIBRARY_PATH:${DIRECTORIES[BUILD]}/lib\""
    full_dotfile_string="\n# --- [Added by MFC | Start Section] --- #\nif [[ -z \"\${MFC_ENV_SH_HEADER_GUARD}\" ]]; then \n\t$export_cmds_0 \n\t$export_cmds_1 \nfi \n# --- [Added by MFC | End Section]   --- #\n"

    i=1
    for dotfile_name in "${DOTFILE_FILENAMES[@]}"; do
        dotfile_path="$HOME/$dotfile_name"

        line_color="${COLORS[GREEN]}"
        status_text="Present"

        if [[ -a "$dotfile_path" ]]; then
            found_dotfile_count=$((found_dotfile_count+1))
            
            if [ $i -ne "1" ]; then
                if [ $i -ne "${#DOTFILE_FILENAMES[@]}" ]; then
                    found_dotfile_list_string="$found_dotfile_list_string, "
                else
                    found_dotfile_list_string="$found_dotfile_list_string, and "
                fi
            fi

            found_dotfile_list_string="$found_dotfile_list_string$dotfile_path"
        else
            line_color=${COLORS[ORANGE]}
            status_text="Absent"
        fi
        
        printf "\r$line_color| %5s | %-70s | %-7s |${COLORS[NONE]}"         \
               "$i/${#DOTFILE_FILENAMES[@]}" "$dotfile_path" "$status_text"

        sleep $SLEEP_DURATION

        i=$((i+1))
    done

    clear_line

    if [ "$found_dotfile_count" -eq "0" ]; then
        echo -e "${COLORS[RED]}"
        print_line_of_char "="
        print_bounded_line "[ERROR] Could not find any dotfiles where we could export the path to the installed libraries."
        print_line_of_char "="
        echo -e "${COLORS[NONE]}"
        exit 1
    fi

    echo -e "\r|--> ${COLORS[GREEN]}Updated $found_dotfile_count dotfiles${COLORS[NONE]}: $found_dotfile_list_string"

else
    echo -e "|--> ${COLORS[GREEN]}The MFC header guard is already present${COLORS[NONE]}: ${COLORS[ORANGE]}No library path will be exported."${COLORS[NONE]}
fi

log_filepath="${DIRECTORIES[LOG]}/FFTW3.log"
cd "${DIRECTORIES[SRC]}/FFTW3"
    show_command_running "|--> (1/$N_DEPENDENCIES) FFTW3: (1/3) Configuring..." \
        log_command "$log_filepath"                                             \
            ./configure --prefix="${DIRECTORIES[BUILD]}"                        \
            --enable-threads                                                    \
            --enable-mpi
    clear_line

    show_command_running "|--> (1/$N_DEPENDENCIES) FFTW3: (2/3) Building..." \
        log_command "$log_filepath"                                          \
            make -j "$JOB_COUNT" CFLAGS="$CFLAGS" CPPFLAGS="$CPPFLAGS"
    clear_line

    show_command_running "|--> (1/$N_DEPENDENCIES) FFTW3: (3/3) Installing..." \
        log_command "$log_filepath"                                            \
            make install
    clear_line

log_filepath="${DIRECTORIES[LOG]}/HDF5.log"
cd "${DIRECTORIES[SRC]}/HDF5"
    show_command_running "|--> (2/$N_DEPENDENCIES) HDF5: (1/3) Configuring..." \
        log_command "$log_filepath"                                            \
            ./configure --enable-parallel                                      \
                        --enable-deprecated-symbols                            \
                        --prefix="${DIRECTORIES[BUILD]}"                       \
                        CC="${USED_UTILITIES[CC]}"                             \
                        CXX="${USED_UTILITIES[CPPC]}"
    clear_line

    show_command_running "|--> (2/$N_DEPENDENCIES) HDF5: (2/3) Building..." \
        log_command "$log_filepath"                                         \
            make -j "$JOB_COUNT" CFLAGS="$CFLAGS" CPPFLAGS="$CPPFLAGS"
    clear_line

    show_command_running "|--> (2/$N_DEPENDENCIES) HDF5: (3/3) Installing..." \
        log_command "$log_filepath"                                           \
            make install prefix="${DIRECTORIES[BUILD]}"
    clear_line

log_filepath="${DIRECTORIES[LOG]}/SILO.log"
cd "${DIRECTORIES[SRC]}/SILO"
    export PYTHON=python3
    export PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $(python3-config --includes) $(python3-config --libs)"

    show_command_running "|--> (3/$N_DEPENDENCIES) SILO: (1/3) Configuring..." \
        log_command "$log_filepath"                                            \
            ./configure --prefix="${DIRECTORIES[BUILD]}"                       \
                        --disable-hzip                                         \
                        --disable-fpzip                                        \
                        --disable-silex                                        \
                        --enable-pythonmodule                                  \
                        --enable-optimization                                  \
                        --with-hdf5="${DIRECTORIES[BUILD]}/include","${DIRECTORIES[BUILD]}/lib" \
                        FC="${USED_UTILITIES[FC90]}"                           \
                        F77="${USED_UTILITIES[FC77]}"                          \
                        CC="${USED_UTILITIES[CC]}"                             \
                        CXX="${USED_UTILITIES[CPPC]}"
    clear_line

    show_command_running "|--> (3/$N_DEPENDENCIES) SILO: (2/3) Building..." \
        log_command "$log_filepath"                                         \
            make -j "$JOB_COUNT" CFLAGS="$CFLAGS" CPPFLAGS="$CPPFLAGS"
    clear_line

    show_command_running "|--> (3/$N_DEPENDENCIES) SILO: (3/3) Installing..." \
        log_command "$log_filepath"                                           \
            make install prefix="${DIRECTORIES[BUILD]}"                       \
    clear_line

echo -e "\r|--> ${COLORS[GREEN]}Successfully built all $N_DEPENDENCIES dependencies${COLORS[NONE]}."

print_logo

if [ "$found_dotfile_count" -ne "0" ]; then
    echo -e "${COLORS[ORANGE]}\n[WARNING] MFC's dependency install script added code to $found_dotfile_count dotfiles ($found_dotfile_list_string) in order to correctly configure your environement variables (such as LD_LIBRARY_PATH). \n${COLORS[NONE]}"
    echo -e "${COLORS[GREEN]}You are now in a new instance of your default shell ($SHELL) and ready to build & run MFC! \n\n${COLORS[NONE]}"
    exec "$SHELL"
fi
