jobid=$(sbatch --parsable <<'SLURM'
#SBATCH --account=b1042
#SBATCH --partition=data-transfer
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --job-name=igv-bam-server
#SBATCH --chdir=/gpfs/projects/b1042/LauberthLab/BenFolder/macs3_results/tracks
#SBATCH --output=/gpfs/projects/b1042/LauberthLab/BenFolder/igv-%j.log

/software/anaconda3/2022.05/bin/python3 - <<'PYTHON'
import http.server
import mimetypes
import os
import re
import secrets
import socket
from urllib.parse import unquote, urlparse

root = os.path.realpath(os.getcwd())
token = secrets.token_urlsafe(24)

class IGVHandler(http.server.BaseHTTPRequestHandler):
    def do_GET(self):
        self.serve_file(head_only=False)

    def do_HEAD(self):
        self.serve_file(head_only=True)

    def do_OPTIONS(self):
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, HEAD, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Range")
        self.end_headers()

    def send_error_response(self, status, message, size=None):
        self.send_response(status)
        self.send_header("Content-Type", "text/plain")
        self.send_header("Access-Control-Allow-Origin", "*")
        if status == 416 and size is not None:
            self.send_header("Content-Range", f"bytes */{size}")
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write((message + "\n").encode())

    def serve_file(self, head_only=False):
        try:
            request_path = unquote(urlparse(self.path).path)
        except Exception:
            return self.send_error_response(400, "Bad URL")

        prefix = "/" + token + "/"
        if not request_path.startswith(prefix):
            return self.send_error_response(403, "Forbidden")

        relative_path = request_path[len(prefix):]
        filename = os.path.realpath(os.path.join(root, relative_path))

        try:
            if os.path.commonpath([root, filename]) != root:
                return self.send_error_response(403, "Forbidden")
        except ValueError:
            return self.send_error_response(403, "Forbidden")

        if not os.path.isfile(filename):
            return self.send_error_response(404, "File not found")

        size = os.path.getsize(filename)
        start = 0
        end = size - 1
        status = 200

        range_header = self.headers.get("Range")
        if range_header:
            match = re.fullmatch(r"bytes=(\d*)-(\d*)", range_header)

            if not match or not any(match.groups()):
                return self.send_error_response(416, "Invalid range", size)

            first, last = match.groups()

            try:
                if not first:
                    suffix_length = int(last)
                    if suffix_length <= 0:
                        raise ValueError
                    start = max(size - suffix_length, 0)
                else:
                    start = int(first)
                    if last:
                        end = min(int(last), size - 1)
            except ValueError:
                return self.send_error_response(416, "Invalid range", size)

            if start < 0 or start >= size or end < start:
                return self.send_error_response(416, "Invalid range", size)

            status = 206

        content_type = (
            mimetypes.guess_type(filename)[0] or "application/octet-stream"
        )

        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(end - start + 1))
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Cache-Control", "no-cache")

        if status == 206:
            self.send_header(
                "Content-Range", f"bytes {start}-{end}/{size}"
            )

        self.end_headers()

        if head_only:
            return

        remaining = end - start + 1

        with open(filename, "rb") as bam_file:
            bam_file.seek(start)

            while remaining:
                chunk = bam_file.read(min(1024 * 1024, remaining))
                if not chunk:
                    break
                self.wfile.write(chunk)
                remaining -= len(chunk)

    def log_message(self, message_format, *args):
        pass

server = http.server.ThreadingHTTPServer(("0.0.0.0", 0), IGVHandler)
port = server.server_address[1]

print(
    f"READY node={socket.gethostname()} port={port} token={token}",
    flush=True
)

server.serve_forever()
PYTHON
SLURM
)

igv_log="/gpfs/projects/b1042/LauberthLab/BenFolder/igv-${jobid}.log"

echo "Submitted job $jobid; waiting for the server..."

while true; do
    igv_ready=$(grep -m1 '^READY ' "$igv_log" 2>/dev/null || true)

    if [[ -n "$igv_ready" ]]; then
        break
    fi

    if ! squeue -h -j "$jobid" | grep -q .; then
        echo "The server job stopped unexpectedly:"
        sed -n '1,160p' "$igv_log" 2>/dev/null
        break
    fi

    sleep 2
done

if [[ -n "$igv_ready" ]]; then
    igv_node=$(sed -n 's/.*node=\([^ ]*\).*/\1/p' <<<"$igv_ready")
    igv_port=$(sed -n 's/.*port=\([^ ]*\).*/\1/p' <<<"$igv_ready")
    igv_token=$(sed -n 's/.*token=\([^ ]*\).*/\1/p' <<<"$igv_ready")

    echo
    echo "JOB ID: $jobid"
    echo
    echo "RUN ON YOUR MAC:"
    echo "ssh -N -o ExitOnForwardFailure=yes -L 60151:${igv_node}:${igv_port} nqp9093@login.quest.northwestern.edu"
    echo
    echo "bam URL:"
    echo "http://127.0.0.1:60151/${igv_token}/H4K16ac_Scrameble_ChIP1_S155_L003_sorted.bw"
    echo
    echo "INDEX URL:"
    echo "http://127.0.0.1:60151/${igv_token}/4K16ac_Scrameble_ChIP1_S155_L003_sorted.bw"
fi


