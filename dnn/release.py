import tarfile
import os
_VERSION_ = '0.3.0.0'
_COMMIT_ = 'na'
tar = tarfile.open("../ochem_dnn_v_{}_{}.tar.gz".format(_VERSION_,_COMMIT_), "w:gz")
__FILES__ = ['apply.py','evaluators.py','models.py','providers.py','train.py','release.py','robust.py','distributed.py']
[tar.add(f) for f in __FILES__]
tar.close()
with tarfile.open("../ochem_dnn_v_{}_{}_data4test.tar.gz".format(_VERSION_,_COMMIT_), "w:gz") as tar:
    [tar.add(source_dir, arcname=os.path.basename(source_dir)) for source_dir in ['./class_only','./reg_only','./mixed']]
# tar = tarfile.open("../ochem_dnn_v_{}_{}_data4testing.tar.gz".format(_VERSION_,_COMMIT_), "w:gz")
# __FILES__ = ['class_only/*.libsvm','reg_only/*.libsvm','mixed/*.libsvm']
# [tar.add(f) for f in __FILES__]
# tar.close()


