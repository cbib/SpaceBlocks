import glob
import os
import shutil


GEOJSON_DIR = config["geojson_path"]
POST_PROCESSING_OUTDIR = config["post_processing_outdir"]


rule copy_geojson:
    """
    Copy the user-provided QuPath GeoJSON for each sample into the workflow output.

    This is a branch after validate_input:
        validate_input -> copy_geojson

    Missing or ambiguous GeoJSON files are reported but do not stop the workflow.
    """
    input:
        validation=rules.validate_input.output.report,

    output:
        done=f"{POST_PROCESSING_OUTDIR}/geojson/{{sample}}.geojson_copied",

    log:
        f"{POST_PROCESSING_OUTDIR}/logs/{{sample}}/copy_geojson.log",

    params:
        geojson_dir=GEOJSON_DIR,
        output_dir=f"{POST_PROCESSING_OUTDIR}/geojson",

    run:
        matches = glob.glob(
            os.path.join(
                params.geojson_dir,
                f"{wildcards.sample}*.geojson"
            )
        )

        os.makedirs(params.output_dir, exist_ok=True)

        with open(log[0], "w") as logfile:

            if len(matches) == 0:
                logfile.write(
                    f"ERROR: No GeoJSON found for sample "
                    f"{wildcards.sample} in {params.geojson_dir}\n"
                )

            elif len(matches) > 1:
                logfile.write(
                    f"ERROR: Multiple GeoJSON files found for sample "
                    f"{wildcards.sample}:\n"
                )
                for path in matches:
                    logfile.write(f"  {path}\n")

            else:
                source = matches[0]
                destination = os.path.join(
                    params.output_dir,
                    os.path.basename(source),
                )

                shutil.copy2(source, destination)

                logfile.write(
                    f"Copied:\n"
                    f"  {source}\n"
                    f"to:\n"
                    f"  {destination}\n"
                )

        # Always create the marker so missing/ambiguous GeoJSON
        # does not block the workflow.
        open(output.done, "w").close()