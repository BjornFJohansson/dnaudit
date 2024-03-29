import os
import uuid
import urllib
import webapp2
import StringIO
import zipfile
import cloudstorage as gcs

test = False

if test:
    class PathWay():
        def __init__(self, sequence_file):
            self.sequence_file = sequence_file
        def zs(self):    
            zipstream = StringIO.StringIO()            
            with zipfile.ZipFile(zipstream, mode='w') as myzip:
                myzip.writestr("hej.txt", self.sequence_file)            
            zipstream.seek(0)
            return zipstream.getvalue()
else:
    from ypkpathway import PathWay

my_default_retry_params = gcs.RetryParams(initial_delay=0.2,
                                          max_delay=5.0,
                                          backoff_factor=2,
                                          max_retry_period=15)
gcs.set_default_retry_params(my_default_retry_params)

from google.appengine.runtime    import DeadlineExceededError
from google.appengine.ext        import blobstore
from google.appengine.ext.webapp import blobstore_handlers
from google.appengine.ext.webapp import template

class MainHandler(webapp2.RequestHandler):
  def get(self):   
    template_file = os.path.join(os.path.dirname(__file__), 'index.html') # <input type="file" name="sequence_file">
    vars_ = { 'upload_url': blobstore.create_upload_url('/upload') }
    rendered = template.render(template_file, vars_)    
    self.response.out.write(unicode(rendered))

class UploadHandler(blobstore_handlers.BlobstoreUploadHandler):
  def post(self):
    try:
        sequence_file   = self.get_uploads('sequence_file')[0].open().read()    
    except IndexError:
        return
               
    write_retry_params = gcs.RetryParams(backoff_factor=1.1)
                                                 
    bucket   = "bukkket"
    filename = "/" + bucket + '/' + str(uuid.uuid4())
   
    gcs_file = gcs.open( filename,
                         'w',
                         content_type='application/zip',
                         options={'x-goog-meta-foo': 'foo',
                                  'x-goog-meta-bar': 'bar'},
                         retry_params=write_retry_params )
    try:                     
        blob_key = blobstore.create_gs_key('/gs' + filename)                     
        result = PathWay(sequence_file).zs()                    
        gcs_file.write(result)
        gcs_file.close()    
        self.redirect('/download/{}'.format(blob_key))
    except DeadlineExceededError:
        self.redirect('/timeout')
            
class TimeOutHandler(webapp2.RequestHandler):
    def get(self):            
        self.response.clear()
        self.response.set_status(500)
        self.response.out.write("This operation could not be completed in time...")    
        
class DownloadHandler(webapp2.RequestHandler):
    def get(self, blob_key):
        template_file = os.path.join(os.path.dirname(__file__), 'download.html')
        vars_ = { 'url': blob_key }
        rendered = template.render(template_file, vars_)
        self.response.out.write(unicode(rendered))

class ServeHandler(blobstore_handlers.BlobstoreDownloadHandler):
    def get(self, resource):
        self.send_blob(resource, save_as="ypk_assembly.zip")        

app = webapp2.WSGIApplication([('/', MainHandler),
                               ('/upload', UploadHandler),
                               ('/download/([^/]+)?', DownloadHandler),
                               ('/serve/([^/]+)?'   , ServeHandler),
                               ('/timeout', TimeOutHandler)],
                               debug=True)
