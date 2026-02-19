#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
while [[ "$root" != "/" && ! -f "$root/CMakeLists.txt" ]]; do
    root="$(dirname "$root")"
done
if [[ ! -f "$root/CMakeLists.txt" ]]; then
    echo "cli.sh: can't find project root (CMakeLists.txt)" >&2
    exit 2
fi

bdir="$root/build"
adc="$bdir/phazeride_adc"
tick="$bdir/parametron_tick"

need_build() {
    [[ -x "$adc" && -x "$tick" ]]
}

do_build() {
    mkdir -p "$bdir"
    cmake -S "$root" -B "$bdir" >/dev/null
    cmake --build "$bdir" -j >/dev/null
}

# defaults
width=8
model="tick"        # tick|logic
signed_mode="auto"  # auto|0|1  (affects dec printing; V computed anyway)
fmt="all"           # dec|hex|bin|all
maj_delay=1
not_delay=1
ticks=0             # 0 -> auto
p_flip="0.0"
seed="-1"           # -1 -> random time seed inside programs
compare="off"       # on|off (only for model=tick)
sat="off"           # off|unsigned|signed
json="0"            # 0|1 (only for adc output)

print_help() {
cat <<'H'
Commands:
  help                         show this help
  quit | exit                  leave
  show                         show current settings
  build                        (re)build binaries into ./build

Settings (either "set key val" or "key val"):
  width N                      1..64
  model tick|logic
  fmt dec|hex|bin|all
  signed auto|0|1
  maj_delay N                  tick only
  not_delay N                  tick only
  ticks N                      0=auto
  p_flip X                     tick only (e.g. 0.001)
  seed N                       -1=auto
  compare on|off               tick only (adc: compare against logic)
  sat off|unsigned|signed      adc only
  json 0|1                     adc only (json-only output)

Actions:
  add A B [CIN]                run phazeride_adc (CIN default 0)
  trace A B [CIN]              run parametron_tick --trace 1 (wave + carry chain)
  A B [CIN]                    shorthand = add

Notes:
  - Numbers: decimal, hex (0x..), negative ok.
  - compare only makes sense for model=tick (adc).
H
}

show_cfg() {
    echo "width=$width model=$model fmt=$fmt signed=$signed_mode json=$json sat=$sat compare=$compare"
    echo "tick: maj_delay=$maj_delay not_delay=$not_delay ticks=$ticks p_flip=$p_flip seed=$seed"
}

run_adc() {
    local a="$1" b="$2" cin="${3:-0}"

    if ! need_build; then
        echo "[cli] build missing -> building..." >&2
        do_build
    fi

    local args=(
        --a "$a" --b "$b" --cin "$cin"
        --width "$width" --model "$model"
        --fmt "$fmt"
        --maj_delay "$maj_delay" --not_delay "$not_delay"
        --ticks "$ticks" --p_flip "$p_flip" --seed "$seed"
        --trace 0
    )

    if [[ "$signed_mode" == "auto" ]]; then
        if [[ "$a" == -* || "$b" == -* ]]; then
            args+=(--signed 1)
        else
            args+=(--signed 0)
        fi
    else
        args+=(--signed "$signed_mode")
    fi

    if [[ "$compare" == "on" ]]; then
        args+=(--compare logic)
    fi

    if [[ "$sat" != "off" ]]; then
        args+=(--sat "$sat")
    fi

    if [[ "$json" == "1" ]]; then
        args+=(--json 1)
    fi

    "$adc" "${args[@]}"
}

run_trace() {
    local a="$1" b="$2" cin="${3:-0}"

    if ! need_build; then
        echo "[cli] build missing -> building..." >&2
        do_build
    fi

    "$tick" \
        --a "$a" --b "$b" --cin "$cin" \
        --maj_delay "$maj_delay" --not_delay "$not_delay" \
        --ticks "$ticks" --p_flip "$p_flip" --seed "$seed" \
        --trace 1
}

set_kv() {
    local k="$1" v="$2"
    case "$k" in
        width) width="$v" ;;
        model) model="$v" ;;
        fmt) fmt="$v" ;;
        signed) signed_mode="$v" ;;
        maj_delay) maj_delay="$v" ;;
        not_delay) not_delay="$v" ;;
        ticks) ticks="$v" ;;
        p_flip) p_flip="$v" ;;
        seed) seed="$v" ;;
        compare) compare="$v" ;;
        sat) sat="$v" ;;
        json) json="$v" ;;
        *) echo "unknown key: $k" >&2; return 1 ;;
    esac
}

validate_cfg() {
    # minimal sanity (keep it simple)
    if ! [[ "$width" =~ ^[0-9]+$ ]] || (( width < 1 || width > 64 )); then
        echo "bad width: $width" >&2; return 1
    fi
    if [[ "$model" != "tick" && "$model" != "logic" ]]; then
        echo "bad model: $model" >&2; return 1
    fi
    case "$fmt" in dec|hex|bin|all) ;; *) echo "bad fmt: $fmt" >&2; return 1 ;; esac
    case "$signed_mode" in auto|0|1) ;; *) echo "bad signed: $signed_mode" >&2; return 1 ;; esac
    case "$compare" in on|off) ;; *) echo "bad compare: $compare" >&2; return 1 ;; esac
    case "$sat" in off|unsigned|signed) ;; *) echo "bad sat: $sat" >&2; return 1 ;; esac
    case "$json" in 0|1) ;; *) echo "bad json: $json" >&2; return 1 ;; esac
    return 0
}

echo "phaserider repl. type 'help'."
show_cfg

while true; do
    IFS= read -r -p "> " line || break
    # trim
    line="${line#"${line%%[![:space:]]*}"}"
    line="${line%"${line##*[![:space:]]}"}"
    [[ -z "$line" ]] && continue

    # split
    read -ra a <<<"$line"
    cmd="${a[0]}"

    case "$cmd" in
        help) print_help ;;
        quit|exit) break ;;
        show) show_cfg ;;
        build) do_build; echo "ok" ;;
        set)
            if (( ${#a[@]} != 3 )); then
                echo "usage: set key value" >&2
                continue
            fi
            set_kv "${a[1]}" "${a[2]}" || continue
            validate_cfg || continue
            ;;
        add)
            if (( ${#a[@]} < 3 || ${#a[@]} > 4 )); then
                echo "usage: add A B [CIN]" >&2
                continue
            fi
            run_adc "${a[1]}" "${a[2]}" "${a[3]:-0}"
            ;;
        trace)
            if (( ${#a[@]} < 3 || ${#a[@]} > 4 )); then
                echo "usage: trace A B [CIN]" >&2
                continue
            fi
            run_trace "${a[1]}" "${a[2]}" "${a[3]:-0}"
            ;;
        width|model|fmt|signed|maj_delay|not_delay|ticks|p_flip|seed|compare|sat|json)
            if (( ${#a[@]} != 2 )); then
                echo "usage: $cmd value" >&2
                continue
            fi
            set_kv "$cmd" "${a[1]}" || continue
            validate_cfg || continue
            ;;
        *)
            # shorthand: "A B [CIN]" => add
            if (( ${#a[@]} == 2 || ${#a[@]} == 3 )); then
                run_adc "${a[0]}" "${a[1]}" "${a[2]:-0}"
            else
                echo "unknown command: $cmd (type help)" >&2
            fi
            ;;
    esac
done
