#!/usr/bin/env bash

# Constants

SCRIPT_FILENAME="$1"
SPIN_ANIMATION=("-" "\\" "|" "/")
SPIN_ANIMATION_WAIT=0.25
SEHLL_TEST_COMMANDS=("declare -A __shell_test_dict__")

# Logo generated here:
# https://patorjk.com/software/taag/#p=moreopts&c=echo&f=Isometric2&t=MFC.
MFC_LOGO=$(cat <<-END
\033[0;34m      ___           ___           ___ 
     /\\  \\         /\\__\\         /\\__\\      \033[0mMulti-component Flow Code\033[0;34m 
    |##\\  \\       /#/ _/_       /#/  / 
    |#|#\\  \\     /#/ /\\__\\     /#/  / 
  __|#|\#\\  \\   /#/ /#/  /    /#/  /  ___ 
 /####|_\#\\__\\ /#/_/#/  /    /#/__/  /\\__\\  \033[0;33m./$SCRIPT_FILENAME -h for help\033[0m 
\033[0;34m \\#\\~~\\  \\/__/ \\#\\/#/  /     \\#\\  \\ /#/  /  
  \\#\\  \\        \\##/__/       \\#\\  /#/  / 
   \\#\\  \\        \\#\\  \\        \\#\\/#/  / 
    \\#\\__\\        \\#\\__\\        \\##/  / 
     \\/__/         \\/__/         \\/__/ \033[0m 
END
)

# Helper Functions

function print_centered_text() {
    if [ ! -z ${TERM+x} ]; then
        printf "%*s\n" $(((${#1}+$(tput cols))/2)) "$1"
    else
        echo "$1"
    fi
}

function print_line_of_char() {
    if [ ! -z ${TERM+x} ]; then
        printf "%$(tput cols)s\n" | tr " " "$1"
    else
        echo "$1$1$1$1$1$1$1$1$1$1$1$1$1$1$1$1"
    fi
}

function print_bounded_line() {
    if [ ! -z ${TERM+x} ]; then
        printf "|| %-$(($(tput cols)-6))s ||" "$1"
    else
        echo "|| $1 ||"
    fi
}

function clear_line() {
    echo -en "\r\033[2K"
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
            sleep $SPIN_ANIMATION_WAIT
        done
    done

    # Once the command has finished running, we run `wait` to get its exit code
    wait $pid
}

function log_command() {
    log_filepath="$1" ; shift
    "$@" >> "$log_filepath" 2>&1
}

# Print Default Text

echo -e "$MFC_LOGO\n"

# Check Shell Compatibility

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
echo -e "\r|--> \033[0;32mSuccessfully ran shell compatibility checks.\033[0m\n"

# Error Handling

function on_exit() {
    exit_code=$?

    if [[ exit_code -ne "0" ]]; then
        echo -e "\033[0;31m"
        print_line_of_char "="
        print_bounded_line "Fatal Error"
        print_bounded_line "---> Exit Code ($exit_code)."
        print_bounded_line "---> Take a look at the output above."
        print_bounded_line "---> Please review the documentation."
        print_line_of_char "="
        echo -e "\033[0m"
        exit $exit_code
    fi
}

trap on_exit EXIT ; set -e ; set -o pipefail ; set -o errtrace

# Define Colors

declare -A COLORS=( [RED]="\033[0;31m"    [GREEN]="\033[0;32m" \
                    [ORANGE]="\033[0;33m" [NONE]="\033[0m"     )
