import subprocess
import threading
import queue
import os

class Worker(threading.Thread):
  def __init__(self, compost_name, jobs, fails, *args, **kwargs):
    self.compost_name = compost_name
    self.jobs = jobs
    self.fails = fails
    super().__init__(*args, **kwargs)
  def run(self):
    while True:
      try:
        work = self.jobs.get(timeout=3)  # 3s timeout
        retorno = export_ligante_to_pdb(self.compost_name, work)
        if (not os.path.exists(self.compost_name)):
          os.makedirs(self.compost_name)
        export_ligante_to_gauss_input(self.compost_name, work)
        if (retorno.returncode != 0):
          self.fails.put({work[0]: work[1]})
      except queue.Empty:
        return
      self.jobs.task_done()

def export_ligante_to_pdb(compost_name, work): #idPai, idProc
  if (not os.path.exists(f"{compost_name}/{work[0]}")):
    os.makedirs(f"{compost_name}/{work[0]}")
  instrucao = f"obabel {work[1]} -o pdb -O  {compost_name}/{work[0]}/{work[0]}.pdb --ff GAFF --gen3d -h --minimize"
  retorno = subprocess.run(instrucao, shell=True, capture_output=True)
  return retorno

def export_ligante_to_gauss_input(compost_name, work): #idPai, idProc
  instrucao = f"obabel {work[1]} -o com -O  {compost_name}/{work[0]}/{work[0]}.com --ff GAFF --gen3d -h --minimize"
  subprocess.run(instrucao, shell=True)
  arquivo_gaus_input = open(f"{compost_name}/{work[0]}/{work[0]}.com", "r")
  conteudo_gaus = []
  for linhas in arquivo_gaus_input.readlines():
    conteudo_gaus.append(linhas)
  conteudo_gaus.pop(0)
  conteudo_gaus.pop(0)
  conteudo_gaus.pop(0)
  conteudo_gaus.pop(0)
  arquivo_gaus_input.close()
  arquivo_gaus_input = open(f"{compost_name}/{work[0]}/{work[0]}.com", "w")
  arquivo_gaus_input.write("%nprocshared=8\n")
  arquivo_gaus_input.write("%mem=1024Mb\n")
  arquivo_gaus_input.write((f"%chk={work[0].lower()}.chk \n"))
  arquivo_gaus_input.write("# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb \n")
  arquivo_gaus_input.write("\n")
  arquivo_gaus_input.write(f"{work[0].lower()}")
  arquivo_gaus_input.write("\n")
  for linhas in conteudo_gaus:
    arquivo_gaus_input.write(linhas)
  arquivo_gaus_input.close()


def counting_threads(compost_name, threads_num, jobs, fails):
    for _ in range(threads_num):
        Worker(compost_name, jobs, fails).start()

def generate(compost_name, dataframe, threads):
  jobs = queue.Queue()
  fails = queue.Queue()
  for id, smile in zip(dataframe.loc[:, "chembl_id"], dataframe.loc[:, "canonical_smiles"]):
      jobs.put_nowait([id, "-:'" + smile + "'"])
  counting_threads(compost_name, threads, jobs, fails)
  print(fails)