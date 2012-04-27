# Argument 1: Download URL (SSL is not working)
# Argument 2: Destination path
# Argument 3: MD5-Hash (Optional)
download_file() {
	if [ "$#" -eq "2" ]; then
		echo "file (DOWNLOAD \"$1\" \"$2\" SHOW_PROGRESS)" > download_file.cmake
	elif [ "$#" -eq "3" ]; then
		echo "file (DOWNLOAD \"$1\" \"$2\" EXPECTED_MD5 \"$3\" SHOW_PROGRESS)" > download_file.cmake
	fi
	echo "Downloading file $1 ..."
	cmake -P download_file.cmake
	rm download_file.cmake
}