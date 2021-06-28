import logging
import os
from pprint import pformat

from Bio import SeqIO





def run_ouango_ContigFilter(config, assemblyUtil, kbase_report):
    """
    This example function accepts any number of parameters and returns results in a KBaseReport
    :param params: instance of mapping from String to unspecified object
    :returns: instance of type "ReportResults" -> structure: parameter
    "report_name" of String, parameter "report_ref" of String
    """
    params = config["params"]  
    shared_folder = config["shared_folder"]
    logging.info('Starting run_ouango_ContigFilter function. Params=' + pformat(params))
   
    # Step 1 - Parse/examine the parameters and catch any errors
    # It is important to check that parameters exist and are defined, and that nice error
    # messages are returned to users.  Parameter values go through basic validation when
    # defined in a Narrative App, but advanced users or other SDK developers can call
    # this function directly, so validation is still important.
    logging.info('Validating parameters.')

        
    for name in ['min_length', 'max_length', 'assembly_input_ref', 'workspace_name']:
        if name not in params:
            raise ValueError('Parameter "' + name + '" is required but missing')
    if not isinstance(params['assembly_input_ref'], str) or not len(params['assembly_input_ref']):
        raise ValueError('Pass in a valid assembly reference string')
    if not isinstance(params['min_length'], int) or (params['min_length'] < 0):
        raise ValueError('Min length must be a non-negative integer') 
    if not isinstance(params['max_length'], int) or (params['max_length'] < 0):
        raise ValueError('Max length must be a non-negative integer') 
    output = {}       
    

    # Step 2 - Download the input data as a Fasta and
    # We can use the AssemblyUtils module to download a FASTA file from our Assembly data object.
    # The return object gives us the path to the file that was created.
    logging.info('Downloading Assembly data as a Fasta file.')
    fasta_file = assemblyUtil.get_assembly_as_fasta({'ref': params['assembly_input_ref']})
    

    # Step 3 - Actually perform the filter operation, saving the good contigs to a new fasta file.
    # We can use BioPython to parse the Fasta file and build and save the output to a file.

    parsed_assembly = SeqIO.parse(fasta_file['path'], 'fasta')
    min_length = params['min_length']
    max_length = params['max_length']

    # Keep a list of contigs greater than min_length
    good_contigs = []
    # total contigs regardless of length
    n_total = 0
    # total contigs over the min_length
    n_remaining = 0
    for record in parsed_assembly:
        n_total += 1
        if len(record.seq) >= min_length and len(record.seq) <= max_length:
            good_contigs.append(record)
            n_remaining += 1
    output = {
        'n_total': n_total,
        'n_remaining': n_remaining
    }

    logging.info('Filtered Assembly to ' + str(n_remaining) + ' contigs out of ' + str(n_total))
    filtered_fasta_file = os.path.join(shared_folder, 'filtered.fasta')
    SeqIO.write(good_contigs, filtered_fasta_file, 'fasta')

    workspace_name = params['workspace_name']
    filtered_path = os.path.join(shared_folder, 'filtered.fasta')
    SeqIO.write(good_contigs, filtered_path, 'fasta')
    # Upload the filtered data to the workspace
    new_ref = assemblyUtil.save_assembly_from_fasta({
        'file': {'path': filtered_path},
        'workspace_name': workspace_name,
        'assembly_name': fasta_file['assembly_name']
    })
    output = {
        'n_total': n_total,
        'n_remaining': n_remaining,
        'filtered_assembly_input_ref': new_ref
    }

    text_message = "".join([
    'Filtered assembly to ',
    str(n_remaining),
    ' contigs out of ',
    str(n_total)
    ])
    # Data for creating the report, referencing the assembly we uploaded
    report_data = {
        'objects_created': [
            {'ref': new_ref, 'description': 'Filtered contigs'}
        ],
        'text_message': text_message
    }
    
    # Initialize the report
    
    report = kbase_report.create({
        'report': report_data,
        'workspace_name': workspace_name
    })
    # Return the report reference and name in our results
    output = {
        'report_ref': report['ref'],
        'report_name': report['name'],
        'n_total': n_total,
        'n_remaining': n_remaining,
        'filtered_assembly_input_ref': new_ref
    }

    # Step 4 - Save the new Assembly back to the system
    logging.info('Uploading filtered Assembly data.')
    new_assembly = assemblyUtil.save_assembly_from_fasta({'file': {'path': filtered_fasta_file},
                                                            'workspace_name': workspace_name,
                                                            'assembly_name': fasta_file['assembly_name']
                                                            })


    # Step 5 - Build a Report and return
    reportObj = {
        'objects_created': [{'ref': new_assembly, 'description': 'Filtered contigs'}],
        'text_message': 'Filtered Assembly to ' + str(n_remaining) + ' contigs out of ' + str(n_total)
    }
    report_info = kbase_report.create({'report': reportObj, 'workspace_name': params['workspace_name']})


    # STEP 6: contruct the output to send back
    output = {'report_name': report_info['name'],
                'report_ref': report_info['ref'],
                'assembly_output': new_assembly,
                'n_initial_contigs': n_total,
                'n_contigs_removed': n_total - n_remaining,
                'n_contigs_remaining': n_remaining
                }
    logging.info('returning:' + pformat(output))
    return output
    