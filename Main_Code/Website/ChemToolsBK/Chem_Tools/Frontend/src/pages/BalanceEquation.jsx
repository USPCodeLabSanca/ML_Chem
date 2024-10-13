import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import '../App.css';

const BalanceEquation = () => {
  const [equation, setEquation] = useState('');
  const [balancedEquation, setBalancedEquation] = useState('');

  const handleSubmit = (e) => {
    e.preventDefault();
    // Funcionalidade para balancear a equação será implementada aqui
    setBalancedEquation('A equação balanceada aparecerá aqui');
  };

  return (
    <div className="ferramentas-quimica">
      <h1>Balancear Equação Química</h1>
      <p className="description">Ferramenta essencial para balancear equações químicas</p>
      
      <div className="conteudo">
        <div className="button-grid">
          <div className="category-column">
            <h2>Calculadora</h2>
            <form onSubmit={handleSubmit}>
              <div>
                <label htmlFor="equation">Digite a equação química:</label>
                <input
                  type="text"
                  id="equation"
                  value={equation}
                  onChange={(e) => setEquation(e.target.value)}
                  placeholder="ex., H2 + O2 = H2O"
                  className="aba input-equation"
                />
              </div>
              <button type="submit" className="aba">Balancear Equação</button>
            </form>
            {balancedEquation && (
              <div className="result">
                <h3>Equação Balanceada:</h3>
                <p>{balancedEquation}</p>
              </div>
            )}
          </div>
          <div className="category-column">
            <h2>Descrição do Problema</h2>
            <p>
              O balanceamento de equações químicas é um processo fundamental na química
              que garante que a lei da conservação da massa seja respeitada. Esta ferramenta
              ajuda você a balancear automaticamente suas equações químicas, economizando
              tempo e reduzindo erros em seus cálculos estequiométricos.
            </p>
            <p>
              Para usar, simplesmente insira sua equação química não balanceada no campo
              de entrada e clique em "Balancear Equação". A ferramenta irá calcular os
              coeficientes corretos para cada elemento, garantindo que o número de átomos
              de cada elemento seja igual em ambos os lados da equação.
            </p>
          </div>
        </div>
      </div>
      <Link to="/" className="aba">Voltar para a Página Inicial</Link>
    </div>
  );
};

export default BalanceEquation;

