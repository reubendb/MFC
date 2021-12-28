#!/usr/bin/env bash

if [ ! -f ./install.sh ]; then
    echo -e "\033[0;31mError: You must call this script from within the 'dependencies' folder.\n\033[0m"
    exit 1
fi

# Note: We shouldn't use wildcards here
DO_SHOW_HEADER=true
if [[ $@ == *"--no-header"* ]] || [[ $@ == *"-nh"* ]]; then
    DO_SHOW_HEADER=false
fi

source ../misc/shell_base.sh install.sh $DO_SHOW_HEADER

# Start of install.sh

function show_help() {
    echo -e "Usage: ./install.sh ${COLORS[ORANGE]}[optional options]${COLORS[NONE]}"
    echo "Optional options:"
    echo -e "  ${COLORS[ORANGE]}-h    (--help)                    ${COLORS[NONE]}Prints this message and exits."
    echo -e "  ${COLORS[ORANGE]}-nh   (--no-header)               ${COLORS[NONE]}Does not print the usual MFC header before executing the script."
    echo -e "  ${COLORS[ORANGE]}-c    (--clean)                   ${COLORS[NONE]}Cleans."
    echo -e "  ${COLORS[ORANGE]}-b    (--build)               [X] ${COLORS[NONE]}Builds all dependencies defined in ${COLORS[ORANGE]}X${COLORS[NONE]} (Default=all). ${COLORS[ORANGE]}X${COLORS[NONE]} is a list of comma-separated dependency names."
    echo -e "  ${COLORS[ORANGE]}-j    (--jobs)                [N] ${COLORS[NONE]}Allows for ${COLORS[ORANGE]}N${COLORS[NONE]} concurrent jobs."
    echo -e "  ${COLORS[ORANGE]}-cc   (--c-compiler)          [X] ${COLORS[NONE]}Uses ${COLORS[ORANGE]}X${COLORS[NONE]} as the C compiler."
    echo -e "  ${COLORS[ORANGE]}-cppc (--cpp-compiler)        [X] ${COLORS[NONE]}Uses ${COLORS[ORANGE]}X${COLORS[NONE]} as the CPP compiler."
    echo -e "  ${COLORS[ORANGE]}-fc77 (--fortran-compiler-77) [X] ${COLORS[NONE]}Uses ${COLORS[ORANGE]}X${COLORS[NONE]} as the Fortran 77 compiler."
    echo -e "  ${COLORS[ORANGE]}-fc90 (--fortran-compiler-90) [X] ${COLORS[NONE]}Uses ${COLORS[ORANGE]}X${COLORS[NONE]} as the Fortran 90 compiler."
    echo -e "  ${COLORS[ORANGE]}-cf   (--c-flags)             [X] ${COLORS[NONE]}Uses ${COLORS[ORANGE]}X${COLORS[NONE]} as the C compiler's flags."
    echo -e "  ${COLORS[ORANGE]}-cppf (--cpp-flags)           [X] ${COLORS[NONE]}Uses ${COLORS[ORANGE]}X${COLORS[NONE]} as the CPP compiler's flags."
    echo -e "  ${COLORS[ORANGE]}-ff   (--fortran-flags)       [X] ${COLORS[NONE]}Uses ${COLORS[ORANGE]}X${COLORS[NONE]} as the Fortran compiler's flags."

    echo -e "\nExample: ./install.sh -b fftw3,hdf5,silo -j 8 -cc mpicc"

    echo -e "\nProvided courtesy of MFC (https://github.com/MFlowCode/MFC).\n"
}

if was_shell_option_used "$*" "-h"; then
    show_help
    exit 0
fi

# Fetch submodules

show_command_running "|--> Fetching submodules... "    \
                     git submodule update --init --recursive
clear_line
echo -e "\r|--> ${COLORS[GREEN]}Fetched submodules.${COLORS[NONE]}"

# CONSTANTS & GLOBAL VARIABLES

declare -A DOWNLOADED_DEPENDENCIES DEFAULT_UTILITIES \
           CHECK_UTILITIES         DIRECTORIES       USED_UTILITIES

DOWNLOADED_DEPENDENCIES=( ["FFTW3"]="http://www.fftw.org/fftw-3.3.10.tar.gz" )

DEFAULT_UTILITIES=( [GIT]="git"         [MAKE]="make"                     \
                    [WGET]="wget|curl"  [TAR]="tar"                       \
                    [PYTHON3]="python3" [PYTHON3-CONFIG]="python3-config" \
                    [CC]="mpicc"        [CPPC]="mpicxx"                   \
                    [FC77]="mpif77"     [FC90]="mpif90"                   )

DOTFILE_FILENAMES=( ".bashrc" ".bash_profile" ".bash_login" ".profile" \
                    ".zshrc"  ".cshrc"        ".zprofile"              )

DIRECTORIES=( [DEPENDENCIES]="$(pwd)" [SRC]="$(pwd)/src" \
              [BUILD]="$(pwd)/build"  [LOG]="$(pwd)/log" )

JOB_COUNT=1

CFLAGS=""
CPPFLAGS=""
FFLAGS=""

DEPENDENCY_NAMES=()
DEPENDENCY_NAMES_TO_BUILD=()

# Fetch the names of all dependencies (submodules + downloads)

DEPENDENCY_NAMES+=("${!DOWNLOADED_DEPENDENCIES[@]}")

submodule_count="$(git submodule status | wc -l)"
index=0 ; while [[ "$index" != "$submodule_count" ]]; do
    DEPENDENCY_NAMES+=("$(git submodule status | tac | tac | head -n $((index+1)) | tail -n 1 | sed "s/\ /\n/g" | head -n 3 | tail -n 1 | sed "s/src\///g")")

    index=$((index+1))
done

DEPENDENCY_NAMES_TO_BUILD=("${DEPENDENCY_NAMES[@]}")

function clean() {
    directory_paths="${DIRECTORIES[LOG]} ${DIRECTORIES[BUILD]}"

    echo -en "\r|--> Removing $directory_paths..."
    rm -rf $directory_paths

    i=1
    for dependency_name in "${!DOWNLOADED_DEPENDENCIES[@]}"; do
        directory_path="${DIRECTORIES[SRC]}:?/$dependency_name"

        clear_line
        echo -en "\r|--> ($i/${#DOWNLOADED_DEPENDENCIES[@]}) Remvoing $directory_path..."
        
        rm -rf "$directory_path"
        i=$((i+1))
    done

    for dependency_name in "${DEPENDENCY_NAMES[@]}"; do
        if [ -d "src/$dependency_name" ]; then
            clear_line
            echo -en "\r|--> Cleaning src/$dependency_name..."
            
            cd "src/$dependency_name"
            make clean 2>&1 > /dev/null
            cd "../../"
        fi
    done

    clear_line
    echo -e "\r|--> ${COLORS[GREEN]}Cleaned.${COLORS[NONE]}"
}

while [[ $# -gt 0 ]]; do
    option="$1"
    case $option in
        -cc  |--c-compiler         ) CHECK_UTILITIES[CC]="$2"   ; shift ; shift ;;
        -fc77|--fortran-compiler-77) CHECK_UTILITIES[FC77]="$2" ; shift ; shift ;;
        -fc90|--fortran-compiler-90) CHECK_UTILITIES[FC90]="$2" ; shift ; shift ;;
        -cppc|--cpp-compiler       ) CHECK_UTILITIES[CPPC]="$2" ; shift ; shift ;;
        -j   |--jobs               ) JOB_COUNT="$2"             ; shift ; shift ;;
        -cf  |--c-flags            ) CFLAGS="$2"                ; shift ; shift ;;
        -cppf|--cpp-flags          ) CPPFLAGS="$2"              ; shift ; shift ;;
        -ff  |--fortran-flags      ) FFLAGS="$2"                ; shift ; shift ;;
        -nh  |--no-header          )                                      shift ;;
        -c   |--clean              ) clean                              ; shift ;;
        -b   |--build              ) 
            DEPENDENCY_NAMES_TO_BUILD=()

            tmp="$(echo "$2" | sed "s/,/\n/g")"

            i=1
            while [[ "$i" != "$(($(echo "$tmp" | wc -l)+1))" ]]; do
                input_name="$(echo "$tmp" | tac | tac | head -n $i | tail -n 1 | sed "s/[0-9]//g" | tr "[:lower:]" "[:upper:]")"
                found_dependency_name=""

                # Check for a match
                bFoundMatch=false
                for dependency_name in "${DEPENDENCY_NAMES[@]}"; do
                    if [[ "$input_name" == "$(echo "$dependency_name" | sed "s/[0-9]//g" | tr "[:lower:]" "[:upper:]")" ]]; then
                        bFoundMatch=true
                        found_dependency_name="$dependency_name"
                        break
                    fi
                done

                if [ $bFoundMatch = false ]; then
                    echo -e "|--> ${COLORS[RED]}Failed to find dependency \"$dependency_name\".${COLORS[NONE]}"
                    echo -e "|--> ${COLORS[ORANGE]}Please check your <-b|--build> options.${COLORS[NONE]}"
                    exit 1
                fi

                # Check for duplicates
                if [[ ! "${DEPENDENCY_NAMES_TO_BUILD[*]}" =~ "${found_dependency_name}" ]]; then
                    DEPENDENCY_NAMES_TO_BUILD+=("$found_dependency_name")
                fi

                i=$((i+1))
            done

            shift ; shift ;;
        *                          )
            echo -e "|--> ${COLORS[RED]}Error: Unknown command line option \"$option\".${COLORS[NONE]}"
            echo -e "|--> ${COLORS[ORANGE]}Please run ./install.sh <-h|--help> for more information.${COLORS[NONE]}"
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

cd "${DIRECTORIES[SRC]}"
    i=1 ; for dependency_name in "${!DOWNLOADED_DEPENDENCIES[@]}"; do
        dependency_link=${DOWNLOADED_DEPENDENCIES[$dependency_name]}

        base_string="|--> Fetching Archives - ($i/${#DOWNLOADED_DEPENDENCIES[@]}) $dependency_name"

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

echo -en "\r|--> ${COLORS[GREEN]}Fetched the source code of all $N_DEPENDENCIES dependencies (downloads and submodules)${COLORS[NONE]}: "

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

dependency_build_counter=0

if [[ "${DEPENDENCY_NAMES_TO_BUILD[*]}" =~ "FFTW3" ]]; then
    dependency_build_counter=$((dependency_build_counter+1))

    log_filepath="${DIRECTORIES[LOG]}/FFTW3.log"
    
    cd "${DIRECTORIES[SRC]}/FFTW3"
    
    show_command_running "|--> ($dependency_build_counter/${#DEPENDENCY_NAMES_TO_BUILD[@]}) FFTW3: (1/3) Configuring..." \
        log_command "$log_filepath"                                                                                      \
            ./configure --prefix="${DIRECTORIES[BUILD]}"                                                                 \
                        --enable-threads                                                                                 \
                        --enable-mpi
    clear_line

    show_command_running "|--> ($dependency_build_counter/${#DEPENDENCY_NAMES_TO_BUILD[@]}) FFTW3: (2/3) Building..." \
        log_command "$log_filepath"                                                                                   \
            make -j "$JOB_COUNT"                                                                                      \
                FFLAGS="$FFLAGS"                                                                                      \
                CFLAGS="$CFLAGS"                                                                                      \
                CPPFLAGS="$CPPFLAGS"
    clear_line

    show_command_running "|--> ($dependency_build_counter/${#DEPENDENCY_NAMES_TO_BUILD[@]}) FFTW3: (3/3) Installing..." \
        log_command "$log_filepath"                                                                                     \
            make install
    clear_line
fi

if [[ "${DEPENDENCY_NAMES_TO_BUILD[*]}" =~ "HDF5" ]]; then
    dependency_build_counter=$((dependency_build_counter+1))

    log_filepath="${DIRECTORIES[LOG]}/HDF5.log"
    
    cd "${DIRECTORIES[SRC]}/HDF5"
    
    show_command_running "|--> ($dependency_build_counter/${#DEPENDENCY_NAMES_TO_BUILD[@]}) HDF5: (1/3) Configuring..." \
        log_command "$log_filepath"                                                                                     \
            ./configure --enable-parallel                                                                               \
                        --enable-deprecated-symbols                                                                     \
                        --prefix="${DIRECTORIES[BUILD]}"                                                                \
                        CC="${USED_UTILITIES[CC]}"                                                                      \
                        CXX="${USED_UTILITIES[CPPC]}"
    clear_line

    show_command_running "|--> ($dependency_build_counter/${#DEPENDENCY_NAMES_TO_BUILD[@]}) HDF5: (2/3) Building..." \
        log_command "$log_filepath"                                                                                  \
            make -j "$JOB_COUNT"                                                                                     \
                 FFLAGS="$FFLAGS"                                                                                    \
                 CFLAGS="$CFLAGS"                                                                                    \
                 CPPFLAGS="$CPPFLAGS"
    clear_line

    show_command_running "|--> ($dependency_build_counter/${#DEPENDENCY_NAMES_TO_BUILD[@]}) HDF5: (3/3) Installing..." \
        log_command "$log_filepath"                                                                                    \
            make install prefix="${DIRECTORIES[BUILD]}"
    clear_line
fi

if [[ "${DEPENDENCY_NAMES_TO_BUILD[*]}" =~ "SILO" ]]; then
    dependency_build_counter=$((dependency_build_counter+1))

    log_filepath="${DIRECTORIES[LOG]}/SILO.log"
    
    cd "${DIRECTORIES[SRC]}/SILO"
    
    export PYTHON=python3
    export PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $(python3-config --includes) $(python3-config --libs)"

    show_command_running "|--> ($dependency_build_counter/${#DEPENDENCY_NAMES_TO_BUILD[@]}) SILO: (1/3) Configuring..." \
        log_command "$log_filepath"                                                                                     \
            ./configure --prefix="${DIRECTORIES[BUILD]}"                                                                \
                        --enable-pythonmodule                                                                           \
                        --enable-optimization                                                                           \
                        --disable-hzip                                                                                  \
                        --disable-fpzip                                                                                 \
                        FC="${USED_UTILITIES[FC90]}"                                                                    \
                        F77="${USED_UTILITIES[FC77]}"                                                                   \
                        CC="${USED_UTILITIES[CC]}"                                                                      \
                        CXX="${USED_UTILITIES[CPPC]}"                                                                   \
                        --with-hdf5="${DIRECTORIES[BUILD]}/include","${DIRECTORIES[BUILD]}/lib"                         \
                        --disable-silex
    clear_line

    show_command_running "|--> ($dependency_build_counter/${#DEPENDENCY_NAMES_TO_BUILD[@]}) SILO: (2/3) Building..." \
        log_command "$log_filepath"                                                                                  \
            make -j "$JOB_COUNT"
    clear_line

    show_command_running "|--> ($dependency_build_counter/${#DEPENDENCY_NAMES_TO_BUILD[@]}) SILO: (3/3) Installing..." \
        log_command "$log_filepath"                                                                                    \
            make install prefix="${DIRECTORIES[BUILD]}"                                                                \
    clear_line
fi

echo -e "\r|--> ${COLORS[GREEN]}Built all ${#DEPENDENCY_NAMES_TO_BUILD[@]} dependencies${COLORS[NONE]}."

echo -e "$MFC_HEADER"

if [ "$found_dotfile_count" -ne "0" ]; then
    echo -e "${COLORS[ORANGE]}\n[WARNING] MFC's dependency install script added code to $found_dotfile_count dotfiles ($found_dotfile_list_string) in order to correctly configure your environement variables (such as LD_LIBRARY_PATH). \n${COLORS[NONE]}"
    echo -e "${COLORS[GREEN]}You are now in a new instance of your default shell ($SHELL) and ready to build & run MFC! \n${COLORS[NONE]}"
    exec "$SHELL"
fi
