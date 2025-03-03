# function for starting a collapsible section
function start_details () {
  local identifier="${1}"
  local title="${2:-$section_title}"

  echo -e "section_start:`date +%s`:${identifier}[collapsed=true]\r\e[0K${title}"
}

# Function for ending the collapsible section
function end_details () {
  local identifier="${1}"

  echo -e "section_end:`date +%s`:${identifier}\r\e[0K"
}
