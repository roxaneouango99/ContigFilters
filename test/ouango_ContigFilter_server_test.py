# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from ouango_ContigFilter.ouango_ContigFilterImpl import ouango_ContigFilter
from ouango_ContigFilter.ouango_ContigFilterServer import MethodContext
from ouango_ContigFilter.authclient import KBaseAuth as _KBaseAuth

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.WorkspaceClient import Workspace 


class ouango_ContigFilterTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('ouango_ContigFilter'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'ouango_ContigFilter',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = ouango_ContigFilter(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  
        cls.prepareTestData()

    @classmethod
    def prepareTestData(cls):
        """This function creates an assembly object for testing"""
        fasta_content = '>seq1 something soemthing asdf\n' \
                        'agcttttcat\n' \
                        '>seq2\n' \
                        'agctt\n' \
                        '>seq3\n' \
                        'agcttttcatgg'

        filename = os.path.join(cls.scratch, 'test1.fasta')
        with open(filename, 'w') as f:
            f.write(fasta_content)
        assemblyUtil = AssemblyUtil(cls.callback_url)
        cls.assembly_ref = assemblyUtil.save_assembly_from_fasta({
            'file': {'path': filename},
            'workspace_name': cls.wsName,
            'assembly_name': 'TestAssembly'
        })

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    
    def test_run_ouango_ContigFilter_max(self):
            ref = "79/16/1"
            result = self.serviceImpl.run_ouango_ContigFilter(self.ctx, {
                'workspace_name': self.wsName,
                'assembly_input_ref': ref,
                'min_length': 100,
                'max_length': 1000000
            })
            
           
    def my_test_run_ouango_ContigFilter_ok(self): 
        # call your implementation
        ret = self.serviceImpl.run_ouango_ContigFilter(self.ctx, 
                                                {'workspace_name': self.wsName,
                                                 'assembly_input_ref': self.assembly_ref,
                                                 'min_length': 10,
                                                 'max_length' : 100000                                              
                                                })
        # Validate the returned data
        self.assertEqual(ret[0]['n_initial_contigs'], 3)
        self.assertEqual(ret[0]['n_contigs_removed'], 1)
        self.assertEqual(ret[0]['n_contigs_remaining'], 2)
    
    def my_test_run_ouango_ContigFilter_min_len_negative(self): 
        with self.assertRaisesRegex(ValueError, 'Min length must be a non-negative integer'):
            self.serviceImpl.run_ouango_ContigFilter(self.ctx,
                                              {'workspace_name': self.wsName,
                                               'assembly_input_ref': '1/fake/3',
                                               'min_length': '-10',
                                              'max_length': '100000'
                                               })
#The params in min_len did not macth the ones in run_o_c. I had to update them.

    def my_test_run_ouango_ContigFilter_min_len_parse(self):
        with self.assertRaisesRegex(ValueError, 'Min length must be a non-negative integer'):
            self.serviceImpl.run_ouango_ContigFilter(self.ctx,
                                              {'workspace_name': self.wsName,
                                               'assembly_input_ref': '1/fake/3',
                                               'min_length': 'ten',
                                               'max_length': '100000'})
    
            
    
    def test_invalid_params(self):
        impl = self.serviceImpl
        ctx = self.ctx
        ws = self.wsName
        # Missing assembly ref
        with self.assertRaises(ValueError):
            impl.run_ouango_ContigFilter(ctx, {'workspace_name': ws,
                'min_length': 100, 'max_length': 1000000})
        # Missing min length
        with self.assertRaises(ValueError):
            impl.run_ouango_ContigFilter(ctx, {'workspace_name': ws, 'assembly_ref': 'x',
                'max_length': 1000000})
        # Min length is negative
        with self.assertRaises(ValueError):
            impl.run_ouango_ContigFilter(ctx, {'workspace_name': ws, 'assembly_ref': 'x',
                'min_length': -1, 'max_length': 1000000})
        # Min length is wrong type
        with self.assertRaises(ValueError):
            impl.run_ouango_ContigFilter(ctx, {'workspace_name': ws, 'assembly_ref': 'x',
                'min_length': 'x', 'max_length': 1000000})
        # Assembly ref is wrong type
        with self.assertRaises(ValueError):
            impl.run_ouango_ContigFilter(ctx, {'workspace_name': ws, 'assembly_ref': 1,
                'min_length': 1, 'max_length': 1000000})
    
    
    def test_run_ouango_ContigFilter_test_min(self):
        ref = "79/16/1"
        params = {
            'workspace_name': self.wsName,
            'assembly_input_ref': ref,
            'min_length': 200000,
            'max_length': 6000000
        }
        result = self.serviceImpl.run_ouango_ContigFilter(self.ctx, params)
        self.assertEqual(result[0]['n_initial_contigs'], 2)
        self.assertEqual(result[0]['n_contigs_removed'], 1)
        self.assertEqual(result[0]['n_contigs_remaining'], 1)
        self.assertTrue(len(result[0]['report_name']))
        self.assertTrue(len(result[0]['report_ref']))
    
        
    def test_run_ouango_ContigFilter_test_max(self):
        ref = "79/16/1"
        params = {
            'workspace_name': self.wsName,
            'assembly_input_ref': ref,
            'min_length': 100000,
            'max_length': 4000000
        }
        result = self.serviceImpl.run_ouango_ContigFilter(self.ctx, params)
        self.assertEqual(result[0]['n_initial_contigs'], 2)
        self.assertEqual(result[0]['n_contigs_removed'], 1)
        self.assertEqual(result[0]['n_contigs_remaining'], 1)
        self.assertTrue(len(result[0]['report_name']))
        self.assertTrue(len(result[0]['report_ref']))
    
        